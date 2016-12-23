function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = Ambulance(x, runlength, seed, other);
% function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = Ambulance(x, runlength, seed, other);
% x is a vector containing the coordinates of the ambulances, (x1, y1, x2, y2, x3, y3)
% runlength is the number of hours of simulated time to simulate
% seed is the index of the substreams to use (integer >= 1)
% other is not used
% Returns Mean response time, no var or gradient estimates.
%   ****************************************
%   *** Code written by German Gutierrez ***
%   ***         gg92@cornell.edu         ***
%   ***                                  ***
%   *** Updated by Shane Henderson to    ***
%   *** use standard calling and random  ***
%   *** number streams                   ***
%   ****************************************
% Last updated June 11, 2011

nAmbulances=1;                  % # of ambulances

FnGrad = NaN;
FnGradCov = NaN;
constraint = NaN;
ConstraintCov = NaN;
ConstraintGrad = NaN;
ConstraintGradCov = NaN;

if (sum(x < 0) > 0) || (sum(x>1)>0) || (sum(size(x) ~= [1, 2 * nAmbulances])>0) || (runlength <= 0) || (seed <= 0) || (round(seed) ~= seed),
    fprintf('x (row vector with %d components)\nx components should be between 0 and 1\nrunlength should be positive and real,\nseed should be a positive integer\n', nAmbulances*2);
    fn = NaN;
    FnVar = NaN;
else % main simulation

lambda=1/120;                    % rate of call arrival
velfk=60; velsk=40;             % Velocities in kilometers
vf=velfk/30; vs=velsk/30;       % Velocities in terms of the unit square
bases = zeros(nAmbulances, 2);  % Base locations
for i=1:nAmbulances,
    bases(i, :) = x(2*i-1:2*i);
end
mus=45; sigmas=15;              % mean and StdDev. of  service times ~ Gamma.

% Generate new streams for call arrivals, call 
[ArrivalTimeStream, LocStream, SceneTimeStream] = RandStream.create('mrg32k3a', 'NumStreams', 3);

% Set the substream to the "seed"
ArrivalTimeStream.Substream = seed;
LocStream.Substream = seed;
SceneTimeStream.Substream = seed;

% Generate random data
OldStream = RandStream.setGlobalStream(ArrivalTimeStream); % Temporarily store old stream

%Generate vector of call arrival times. (used to calculate service time per call).
InterarrivalTimes = exprnd(1/lambda, 1, runlength);
CallTimes = cumsum(InterarrivalTimes);
Tmax = CallTimes(runlength); % German's code uses the time to finish, in addition to the number of calls

NumCalls = runlength;

% Generate scene times
RandStream.setGlobalStream(SceneTimeStream);
Serv = gamrnd((mus/sigmas)^2,mus/(mus/sigmas)^2, 1, NumCalls);

%Generate Call Locations - use acceptance rejection
RandStream.setGlobalStream(LocStream);
%CallLocations=rand(2,NumCalls);
CallLocations=zeros(2,NumCalls);
i=1;
while i~=NumCalls+1
    u = rand(3,1);
    if 1.6 * u(3) <= 1.6-(abs(u(1)-.8)+abs(u(2)-.8))
        CallLocations(:,i)=u(1:2);
        i=i+1;
    end
end

% Restore old random number stream
RandStream.setGlobalStream(OldStream);

ExitTimes=zeros(1,NumCalls);                 %keeps track of time at which service for each call is completed.
AmbulanceArrivalTimes=zeros(1,NumCalls);     %keeps track of time at which ambulance arrived to each call.


% The following matrix will contain updated information of the location of 
% last call responded by each ambulance, i.e. col.1=Time at which last
% call was finished (travel+service) for ambulance 1, col. 2= X location of 
% last call, col.3= Y location of last call, col. 4= X location of its base 
% and col. 5= Y location of its base.

Ambulances=zeros(nAmbulances,5);
Ambulances(:,4:5)=bases;
Ambulances(:,2:3)=bases;

%Loop through all calls, assign the ambulance that will respond to it and
%update the finish time for the given ambulance. To do this, we must look
%at available ambulances (if any) and assign the closest one. Else, the
%next available ambulance will respond to the call.

%DistanceToCall;     %keeps track of distance between ambulance and present call
%xcurrent;           %current x location of an Ambulance (at time of present call)
%ycurrent;           %current y location of an Ambulance (at time of present call)
minTime=Tmax+10000;  %keeps track of time at which next ambulance will be available(if all currently busy)
%closestA;           %keeps track of closest Ambulance (out of the ones available)
%minDistance;        %minimum distance to call out of available ambulances
%xcall; ycall;       %keeps location of call
%xlc; ylc;           %keeps location of last call serviced
%xb; yb;              %keeps location of an ambulance's base
for i=1:NumCalls
   if(CallTimes(i)~=0)
    closestA=0;
    minDistance=1000000;
    xcall=CallLocations(1,i);
    ycall=CallLocations(2,i);
    for j=1:nAmbulances
        %check if ambulance j is available, if so calculate how far it is
        %from the call. Keep track of closest available ambulance
        if(Ambulances(j,1)<CallTimes(i))
            xlc=Ambulances(j,2);
            ylc=Ambulances(j,3);
            xb=Ambulances(j,4);
            yb=Ambulances(j,5);
            %enough time between current and last call for ambulance to be
            %at base
            if(CallTimes(i)-Ambulances(j,1)>(abs(xlc-xb)+abs(ylc-yb))/vs)
                DistanceToCall= abs(xcall-xb)+abs(ycall-yb);
            %ambulances at y-location of base, still traveling towards
            %x-location of base.
            elseif (CallTimes(i)-Ambulances(j,1)>abs(ylc-yb)/vs)
                %to discover the horizontal direction of movement
                %(right/left).
                    if(xb-xlc>0)
                        xcurrent=xlc+vs*(CallTimes(i)-Ambulances(j,1)-abs(ylc-yb)/vs);
                        DistanceToCall=abs(ycall-yb)+ abs(xcall-xcurrent);
                    else
                        xcurrent=xlc-vs*(CallTimes(i)-Ambulances(j,1)-abs(ylc-yb)/vs);
                        DistanceToCall=abs(ycall-yb)+abs(xcall-xcurrent);
                    end
            %ambulance still trying to get to y-location of base.
            else
                %to discover the vertical direction of movement (up/down)
                    if(yb-ylc>0)
                        ycurrent=ylc+vs*(CallTimes(i)-Ambulances(j,1));
                        DistanceToCall=abs(ycall-ycurrent)+abs(xcall-xlc);
                    else
                        ycurrent=ylc-vs*(CallTimes(i)-Ambulances(j,1));
                        DistanceToCall=abs(ycall-ycurrent)+abs(xcall-xlc);
                    end
            end
            %if ambulance closer than the closest available ambulance so
            %far, keep track of it.
            if(DistanceToCall<minDistance)
                closestA=j;
                minDistance=DistanceToCall;
            end
        end
    end
    %If there is an available ambulance, dispatch it. i.e. set its last
    %call to be this one, update its following finish time, x and y
    %locations of last call serviced.
    if minDistance~=1000000
        ExitTimes(i)=CallTimes(i)+DistanceToCall/vf+Serv(i);
        AmbulanceArrivalTimes(i)= CallTimes(i)+DistanceToCall/vf;
        Ambulances(closestA,1)=CallTimes(i)+DistanceToCall/vf+Serv(i);
        Ambulances(closestA,2)=xcall;
        Ambulances(closestA,3)=ycall;
    else
        %No available ambulances, therefore the next available ambulance 
        %will service the call.
        minTime=Tmax+10000;
        for j=1:nAmbulances
            %find next available ambulance
            if(Ambulances(j,1)< minTime)
                minTime=Ambulances(j,1);
                ClosestA=j;
            end    
        end
        %update the next finish time, x and y locations of last call
        %serviced for ambulance that will respond to the call.
        ExitTimes(i)=minTime+((abs(Ambulances(ClosestA,2)-xcall)+abs(Ambulances(ClosestA,3)-ycall))/vf)+Serv(i);
        AmbulanceArrivalTimes(i)= minTime+((abs(Ambulances(ClosestA,2)-xcall)+abs(Ambulances(ClosestA,3)-ycall))/vf);
        Ambulances(ClosestA,1)=minTime+((abs(Ambulances(ClosestA,2)-xcall)+abs(Ambulances(ClosestA,3)-ycall))/vf)+Serv(i);
        Ambulances(ClosestA,2)=xcall;
        Ambulances(ClosestA,3)=ycall;
    end
    end
end


%Calculate AvgResponseTime and Standard deviation. Use Ambulance arrival
%times to calls - time of call.
fn = mean(AmbulanceArrivalTimes(1,1:nnz(CallTimes)-1)-CallTimes(1,1:nnz(CallTimes)-1));
FnVar = NaN; % Variance estimate of function estimate not available

end % if input parameters are valid


