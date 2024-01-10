function varargout = integrate_with_restarts(odesolver, odefun, tspan, y0, varargin)
% sol = integrate_with_restarts(odesolver, odefun, tspan, y0, varargin)
% [T,Y,Tr,sol] = integrate_with_restarts(odesolver, odefun, tspan, y0, varargin)
% [T,Y,TE,YE,IE,Tr,sol] = integrate_with_restarts(odesolver, odefun, tspan, y0, varargin)
%
% Uses the specified to integrate a function.
% If the integrator is stopped in between (e.g. by an event function),
% the integration is continued at the current point.
%
% INPUT:  odesolver --> handle to ode-solver (e.g. @ode45)
%                       will be called as odesolver(odefun, tspan, y0, options)
%                t0 --> initial time point
%                tf --> final time point
%          varargin --> additional arguments as key-value-pairs as below
%
% Optional input arguments passed as key-value-pairs:
%     solveroptions --> Cell array of arguments passed to the odesolver
%    stateupdatefun --> Function that may update the states at each stop,
%                       and is of form [newstate switchData] = updfun(t,y,sol),
%                       see below for details.
%    modelupdatefun --> Function delivering new right hand sides at each stop
%                       and is of form [newrhs switchData] = modelUpdateFunction(t,y,sol)
%                       see below for details.
%            cadlag --> Flag indicating wether the output shall be made cadlag,
%                       (continue a droite, limite a gauche) by shifting the end
%                       of each interval is shifted by an epsilon to the left.
%                       NOTE: This allows easy usage of "interp" and "deval".
%    diagnosticmode --> Flag to set the diagnostic mode (default: true, false)
%               xE0 --> Use user defined initial event states xE(t=0)
% addSwitchData2sol --> Note that standard optional arguments to the
%                       <odefun> are static - they can not be changed, and 
%                       they are not available in the <modelUpdater>. 
%                       If odefun with (new) optional arguments needs to be 
%                       defined in the <modelUpdater> use the switchData container 
%                       to pass this data - it is available in the sol object. 
%                       The container might also be used to exchange data 
%                       between <modelupdatefun> and <stateupdatefun>, or 
%                       for any other monitoring information
%
% OUTPUT:    T --> vector of times
%            Y --> matrix of results as from odesolver
%           TE --> event times of integator
%           YE --> event states of integrator
%           IE --> event indices of integrator
%           Tr --> time points of integrator restarts
%          sol --> solution structure of integrator additional fields
%                  .restarts  (= Tr = time points of restars)
%                  .intervals (number of integration intervals)
%
%
% State Update Function:
%    The function is of form newstate = stateUpdFun(t,y,sol) where
%         t --> current time
%         y --> current state
%       sol --> current solution structure ('sol' from odesolver or odextend)
%
% Model Update Function:
%    The function is of form newrhs = modelUpdateFunction(t,y,sol)
%    With same arguments as the state update function and delivering
%    a funtion handle newrhs that will be used as right hand side on the
%    next integration interval.
%
%
% (c) Andreas Sommer, Nov2016, Feb2017
% andreas.sommer@iwr.uni-heidelberg.de
% email@andreas-sommer.eu
%
% v10 modified by:  Tilman Barz, 19/08/2018, tilman.barz@ait.ac.at
%               . switchData container added
%               . input analysis added
%               . detailed output during integration added
%               . additional info my_XXX added to output object
%
% v20 modified by:  Tilman Barz, 20/09/2018, tilman.barz@ait.ac.at
%               . detection of duplicated events added
%               . rejection of events and consistency checks added
%
% v22 modified by:  Tilman Barz, 27/09/2018, tilman.barz@ait.ac.at
%               . additional consistency check for detected events added
%
% v23 modified by:  Tilman Barz, 30/01/2019, tilman.barz@ait.ac.at
%               . xE0 output added for initialization of directions 
%
% v23b modified by:  Tilman Barz, 01/08/2019, tilman.barz@ait.ac.at
%               . first version supporting input analyis when 
%                 detecting events crossing with either negativ or
%                 positiv direction only -> needs to be continued...
%                 <CrossingSign> analysis extded
%               . solveroptions now mandatory including handle to event fct
%

  
% NOTES:
%   * If an event function signals a terminal event, matlab *always* adds the
%     terminal time and terminal state to the result vector. This happens also,
%     if that time point was not specifically requested in tspan.





% ===========================
% PROCESSING ARGUMENTS

% set defaults
switchData     = [];
stateupdatefun = [];
modelupdatefun = [];
cadlag = false;
solveroptions  = {};
diagnosticMode = true;  % DEBUG
xE0_option = false;

% deprecated arguments (for backward compatibility
if hasOption(varargin, 'updatefunction')
   warning('Deprecated argument "updatefunction", use "stateupdatefun" instead!')
   stateupdatefun = getOption(varargin, 'updatefunction');
end
if hasOption(varargin, 'interpolatable')
   warning('Deprecated argument "interpolatable", use "cadlag" instead!')
   cadlag = getOption(varargin, 'interpolatable');
end

% process optional arguments
if hasOption(varargin, 'addSwitchData2sol'), switchData     = getOption(varargin, 'addSwitchData2sol'); end
if hasOption(varargin, 'diagnosticmode'),    diagnosticMode = getOption(varargin, 'diagnosticmode'); end
if hasOption(varargin, 'stateupdatefun'),    stateupdatefun = getOption(varargin, 'stateupdatefun'); end
if hasOption(varargin, 'modelupdatefun'),    modelupdatefun = getOption(varargin, 'modelupdatefun'); end
if hasOption(varargin, 'solveroptions'),     solveroptions  = getOption(varargin, 'solveroptions');  
else error('<solveroptions> have to be defined, they need to include the handle to an event function!'); end
if hasOption(varargin, 'cadlag'),            cadlag         = getOption(varargin, 'cadlag');         end
if hasOption(varargin, 'initialEventState'), xE0            = getOption(varargin, 'xE0');            
                                             xE0_option     = true;                                  end

% error checks on number of output arguments:
if nargout==0, warning('No output args specified. Delivering sol structure.'), end
if nargout >8, warning('Too many output args specified! Trying my best...'), end



% END OF PROCESSING ARGUMENTS
% ===========================






% ===========================
% PREPARE INTEGRATION

    % extract start and end time
    t0 = tspan(1);
    tf = tspan(end);

    diagmsg('\n==========================================================================================================\n')
    diagmsg('=> Preparing for integration from time: [%10.4f] -> [%10.4f]\n', t0, tf)
    if cadlag
    diagmsg('\n * Cadlag <on> (right continuous with left limits), allows easy usage of "interp" and "deval" to obtain solution at desired time points. \n');
    end
    
    % prepare integration
    intvl = 1;    % loop counter as Matlab restarts
    my_intvl = 1; % loop counter as real events are found
    ticID = tic;  % start timer
    eventCnt = 0; % how many events have been detected
    rejEvCnt = 0; % how many events have been rejected

    % set initial time and state
    initialTime  = t0;
    initialState = y0;

    % extract info on monitored event states
%     initialEventState   = solveroptions.Events(t0,y0);
    [initialEventState,initialIsterminal,initialDirection] = solveroptions.Events(t0,y0);
    if xE0_option
        if sum(size(xE0)) ~= sum(size(initialEventState))
        diagmsg('<ERROR> User defined <xE0> have different dimension as output from event function. \n');
        error('Stop.')
        end
        initialEventState = xE0;
        diagmsg(' * Using user defined initial event states: xE0. \n');
        diagmsg('   The sign(xE0) is used to initialize the output <my_SgnE> the integration start. \n');
    end

    initialCrossingSign = sign(initialEventState);
    
% check eventStates    
    if find(initialEventState==0)
        tmp1 = find(initialEventState==0);
        diagmsg('<ERROR> The event states with index ['); diagmsg(' %d ',tmp1);
        diagmsg('] are zero!\n')
        diagmsg('This is not a good start/ initialization!  => ')
        error('WARNING')
    end
        diagmsg('   [%3d] events are monitored.\n',length(initialEventState)); 
        diagmsg('      *  Initial event states:                    '); 
        diagmsg('%3d|',initialCrossingSign); diagmsg('\n');
        diagmsg('         (positive sign means "last" crossing was with positive sign/gradient)\n')
% check isTerminal
    if find(initialIsterminal~=1)
        tmp2 = find(initialIsterminal~=1);
        diagmsg('<ERROR> non terminal events detected, event index ['); diagmsg(' %d ',tmp2);
        diagmsg('This is not supported! \n')
        error('stop')
    end
        diagmsg('      *  All events are terminal events.\n'); 
% check direction
        diagmsg('      *  Detection of zero-crossing from:         '); 
        diagmsg('%3d|',initialDirection); diagmsg('\n');
        diagmsg('         (-1 means negative direction only.)\n'); 
        diagmsg('         ( 0 means both directions.)\n'); 
        diagmsg('         (+1 means positive direction only.)\n'); 
    if find(initialDirection~=0)
        diagmsg('<WARNING> direction option <0> only is supported! \n');
%         error('stop')
    end        
    

% ===========================
% FIRST INTEGRATION INTERVAL
    diagmsg('==========================================================================================================\n')
    diagmsg('=> Interval #[%3d] @time: [%10.4f] ->', intvl, initialTime);
sol = odesolver(odefun, [initialTime tf], initialState, solveroptions);
    diagmsg(' [%10.4f]', sol.x(end));

   
   % store infos to be added to sol @t0
    my_xe{my_intvl}   = initialTime;
    my_ye{my_intvl}   = initialState;
    my_ie{my_intvl}   = 0; % no event at t0
    my_SgnE{my_intvl} = initialCrossingSign;
    my_rhsfuncs{my_intvl} = odefun;
    %     restarts = cell(1);
    %     rhsfuncs = cell(1);

    
    % for cadlag only
    restartTimeIndex = length(sol.x);
%     % for loop only
%     currentTime  = sol.x(end);


% ===========================
% INTEGRATION LOOP
while (sol.x(end) < tf) 

    % get infos on triggered events = last stop
    % analyse the last detected events and update 
    diagmsg('     * event detected *\n');
    [eventsFoundInIDX initialCrossingSign eventCnt rejEvCnt] = ...
       analyseLastEventsAndUpdate_my_arrays(sol,eventCnt,rejEvCnt,initialCrossingSign,initialTime);
   
    % Update and initial time and state
    initialTime  = sol.x(end); 
    initialState = sol.y(:,end);
    % checker
    if (sum(abs(sol.y(:,end) - sol.ye(:,end)))> sqrt(eps));   
        error('strange here'); 
    end
    
    % now prepare the start of the next Matlab interval
    intvl = intvl + 1;

    if ~isempty(eventsFoundInIDX) % only if a "valid" event was detected
        
        % now prepare the start of the next Matlab interval
        my_intvl = my_intvl + 1;
        
        % append the info on last detected event(s)
        my_xe{my_intvl}   = sol.xe(end);
        my_ye{my_intvl}   = sol.ye(:,end); % might be updated by < stateupdatefun >
        my_ie{my_intvl}   = eventsFoundInIDX; 
        my_SgnE{my_intvl} = initialCrossingSign;
        my_rhsfuncs{my_intvl} = odefun;    % might be updated by < modelupdatefun >
    
        % STORE additional infos in SOL object
        % they will be made available in the <stateupdatefun> and <modelupdatefun>
        %     sol.restarts = cell2mat(restarts);
        sol.my_xe       = my_xe;
        sol.my_ie       = my_ie; 
        sol.my_SgnE     = my_SgnE;
        sol.my_ye       = my_ye;
        sol.my_rhsfuncs = my_rhsfuncs;
        % -> updates need to be stored in sol, if not they are lost! 
        sol.switchData = switchData;

        % update the initial state for the next interval
        if ~isempty(stateupdatefun)
                diagmsg('   ---> enter the stateupdatefun --->\n');
%                 diagmsg('        to update the state a time point number %3d --->\n',my_intvl);
            [initialState switchData] = stateupdatefun(initialTime, initialState, sol);
                diagmsg('   <-- leave the stateupdatefun <---\n');
            % append the info on last event(s)
            my_ye{my_intvl} = initialState;
            % STORE additional infos in SOL object
            % they will be made available in the <stateupdatefun> and <modelupdatefun>
            sol.my_ye = my_ye;       % might be updated by < stateupdatefun >
            sol.switchData  = switchData;
        else
            % no update...
            sol.my_ye = my_ye; 
        end

        % update the model for the next interval
        if ~isempty(modelupdatefun)
                diagmsg('   ---> enter the modelupdateFun --->\n');
            [odefun switchData] = modelupdatefun(initialTime, initialState, sol);
                diagmsg('   <-- leave the modelupdateFun <---\n');
            % append the info on last event(s)
            my_rhsfuncs{my_intvl} = odefun;
            % STORE additional infos in SOL object
            % they will be made available in the <stateupdatefun> and <modelupdatefun>
            sol.my_rhsfuncs = my_rhsfuncs;
            sol.switchData  = switchData;
        else
            % no update...
            sol.my_rhsfuncs = my_rhsfuncs; 
        end
        
    end
    
    % store initial integration (re-)start time
    %     restarts{intvl} = initialTime;
    %     rhsfuncs{intvl} = odefun;
   
    % extend solution and adjust sol structure
    diagmsg('----------------------------------------------------------------------------------------------------------\n')
    diagmsg(  '=> Interval #[%3d] @time: [%10.4f] ->', intvl, initialTime);
    
    sol = odextend(sol, odefun, tf, initialState, solveroptions);
    diagmsg(  ' [%10.4f]', sol.x(end));

    
    % % %     % Update and store initial time and state
    % % %     initialTime  = sol.x(end);
    % % %     % this has to be done here, because of while loop and cadlag!!!
    
    % possibly adjust solution for being cadlag
    if cadlag
        sol.x(restartTimeIndex) = sol.x(restartTimeIndex) - eps(sol.x(restartTimeIndex));
        % important: needs to be done afterwards, 
        %   sol.x(end)  is monitroed in the while loop
        restartTimeIndex = length(sol.x);
    end
       
end
% stop timer
elapsedTime = toc(ticID);

% ===========================
% PRINT some statistics
diagmsg('\n')
    diagmsg('   stopped @ time: %18.18g <-> t_f: %18.18g \n', sol.x(end),tf);
diagmsg('\n')
diagmsg('   =>=>=>==>>>   F.I.N.I.S.H.E.D   <<<==<=<=<=\n')
diagmsg('   Integration took %g seconds\n', elapsedTime);
diagmsg('   Sum of (initial)start and re-starts: %d\n', intvl);
diagmsg('   Total number of detected events:    %d  (%d have been rejected)\n', eventCnt, rejEvCnt);


    % STORE additional infos in SOL object
    % they will be made available in the <stateupdatefun> and <modelupdatefun>
    %     sol.restarts = cell2mat(restarts);
    sol.my_xe       = my_xe;
    sol.my_ye       = my_ye;       % might be updated by < stateupdatefun >
    sol.my_ie       = my_ie; 
    sol.my_SgnE     = my_SgnE;
    sol.my_rhsfuncs = my_rhsfuncs; % might be updated by < modelupdatefun >
    
    % add additional information to sol
    % might be used in < stateupdatefun > or < modelupdatefun >
    sol.switchData = switchData;


% ===========================
% ADJUSTMENTS

% NOT NEEDED ANYMORE --> ALREADY DONE IN INTEGRATION LOOP
% make interpolatable function if requested (ensures strong monotonic time)
% NOTE: This makes the function cadlag 
%if cadlag
%   idx = ismember(sol.x, restarts(2:end));     % find restart times
%   idx = find(idx);                            % get the linear indices
%   idx = idx(2:2:end);                         % indices of new interval starts
%   sol.x(idx) = sol.x(idx) + eps(sol.x(idx));  % adjust the respective entries
%end



% ===========================
% OUTPUT ASSEMBLY


    % Assemble output arguments as requested by user
    if (nargout==0) % no output args: warn and return sol structure
       warning('No output args specified. Delivering sol structure.')
       varargout{1} = sol;
    end
    if (nargout == 1) % return sol structure
       varargout{1} = sol;
    end
    if (nargout >= 2) % requested: T, Y, [Tr]
       if length(tspan)==2
          varargout{1} = sol.x;   % if onty t0 and tf specified, 
          varargout{2} = sol.y;   % return the steps as determined by solver
       else
          varargout{1} = tspan;   % otherwise return tspan and 
          varargout{2} = deval(sol,tspan); % interpolate at requested times
       end
       varargout{3} = cell2mat(restarts);
    end
    if (nargout >=4 ) % requested: T, Y, TE, YE, IE, [Tr}
       varargout{3} = sol.xe;
       varargout{4} = sol.ye;
       varargout{5} = sol.ie;
       varargout{6} = cell2mat(restarts);
    end

    % Too manx output arguments?
    maxnargout = 6;
    if (nargout > maxnargout)
       warning('Invalid number of output arguments: %d (max: %d)!', nargout, maxnargout);
    end

diagmsg('==========================================================================================================\n\n')


% FINITO
return


% ===========================================================
% ===========================================================



% === HELPERS ===
   function diagmsg(messagestr, varargin)
      if diagnosticMode
         fprintf(messagestr, varargin{:}) %#ok<PRTCAL>
      end
   end


% === HELPERS ===
function [eventsFoundInIDXRed initialCrossingSign eventCnt rejEvCnt] = ...
        analyseLastEventsAndUpdate_my_arrays(sol,eventCnt,rejEvCnt,initialCrossingSign,initialTime)

    % provide info on last stopping events
    eventsPosi       = [eventCnt+1:length(sol.ie)];
    eventsFoundInIDX = sol.ie(eventsPosi);
    
    % how many new events detected
    eventCnt = eventCnt + length(eventsFoundInIDX);
    % how many events in total?
    % totalNumEvents = length(sol.xe);

    % check if the event indices are unique
    [numRepVal1 tmp1 tmp1] = isUniqueIntegerVector(uint8(eventsFoundInIDX));  
    % initialize the rejector of events
    eventRejector = zeros(size(eventsPosi));

    if numRepVal1>0 
        diagmsg('   [%d] events are detected      ->  !!! multiple events for same event state !!!', ...
                   length(eventsFoundInIDX));
        diagmsg('    | total events now: [%d] \n', eventCnt);               
        flagRepVal = 1;
        
        % which events are too close to the start? 
        diffFromStart = abs(sol.xe(eventsPosi) - initialTime);
        evRejIDX = find( diffFromStart < sqrt(eps(sol.xe(eventsPosi)))) ;
        eventRejector(evRejIDX) = 1;
        
        % check if their are still problems with selected reduced event set
        % check if the selected event indices are unique 
        tmpRedEvFoundInIDX = eventsFoundInIDX(find(~eventRejector));
        [numRepVal2 tmp1 tmp1] = isUniqueIntegerVector(uint8(tmpRedEvFoundInIDX));  
        if numRepVal2>0 
            error('This is an unexpected feature (bug?) in Matlabs event detection !')
        end
        
    else % everything ok!
        diagmsg('   [%d] events reported @time: %18.18f | total events now: [%d] \n', ...
                   length(eventsFoundInIDX),sol.x(end),eventCnt);
        flagRepVal = 0;
    end
    
    % check if selected event times of the reduced event set are unique
    remainingEvPosi = eventsPosi(find(~eventRejector));
    remainingEvTimes = sol.xe(remainingEvPosi);   
%     if length(unique(remainingEvTimes))>1 
    % error('This is an unexpected feature (bug?) in Matlabs event detection !')
        
        % which events are too close to the start? 
        diffFromStart = abs(remainingEvTimes - initialTime);
        evRejIDX = find( diffFromStart < sqrt(eps(remainingEvTimes))) ;
        eventRejector(evRejIDX) = 1;

%     end
        
    diagmsg('       change the crossing sign from: ')
    diagmsg('%3d|', initialCrossingSign); diagmsg('\n')
          
    % go through all events
    for io = 1:length(eventsFoundInIDX)
        chIDXev = eventsFoundInIDX(io);
        diffFromStart = abs(sol.xe(eventsPosi(io)) - initialTime);

        diagmsg('    *  y-IDX: %4d,               to: ',chIDXev);  
        
        if eventRejector(io)
            % reject
            diagmsg('%3d|', initialCrossingSign);
            diagmsg('   @time: %18.18f    !!! Rejected !!!\n', sol.xe(eventsPosi(io)))
        else
            % accept
            % now update the switched els in crossing signs
            te = sol.xe(eventsPosi(io));
            [initialCrossingSign] = updCrossingSigns(te,chIDXev,initialCrossingSign,initialDirection);
%             initialCrossingSign(chIDXev) = - initialCrossingSign(chIDXev);
%             diagmsg('%3d|', initialCrossingSign);
%             diagmsg('   @time: %18.18f\n', sol.xe(eventsPosi(io)))
        end
    end        
    
%                     if length(find(sol.x(end)==sol.xe)) ~= eventsFoundInIDX
                        % cadlag makes this test kaputt
                %             error('events stored by sol object erroneous?')
%                     end

    rejEvCnt = rejEvCnt + sum(eventRejector);
      
    % only unique/valid events will be stored
    eventsFoundInIDXRed = eventsFoundInIDX(find(~eventRejector));
end


% === HELPERS ===
function [initialCrossingSign] = updCrossingSigns(te,chIDXev,initialCrossingSign,initialDirection)
    drct = initialDirection(chIDXev);
    if drct==0 
        initialCrossingSign(chIDXev) = - initialCrossingSign(chIDXev);
        diagmsg('%3d|', initialCrossingSign);
        diagmsg('   @time: %18.18f\n', te);
    elseif drct==-1  
        initialCrossingSign(chIDXev) = - 1;
        diagmsg('%3d|', initialCrossingSign);
        diagmsg('   @time: %18.18f\n', te);
        diagmsg('       (Note: This event is triggered for negative direction only.)\n');
    elseif drct== 1 
        initialCrossingSign(chIDXev) =   1;
        diagmsg('%3d|', initialCrossingSign);
        diagmsg('   @time: %18.18f\n', te);
        diagmsg('       (Note: This event is triggered for positive direction only.)\n');
    else
        error('something wrong');
    end
end



% end of main function
end

