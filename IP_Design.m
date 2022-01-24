%% Inverted Pendulum: Modern Control Approach

%% Inverted Pendulum properties

M = 0.5;                    % mass of cart in kg
m = 0.2;                    % mass of pendulum in kg
b = 0.1;                    % coefficient of friction for cart in N/m/sec 
I = 0.006;                  % mass moment of inertia in kg.m^2
g = 9.8;                    % acceleration due to gravity in m/s^2
l = 0.3;                    % length to pendulum centre of mass in m

%% State-Space Modeling

p=I*(M+m)+(M*m*l^2);

A=[0       1              0         0;    % system matrix
   0 (-(I+m*l^2)*b)/p (m^2*g*l^2)/p 0; 
   0       0              0         1; 
   0 (-m*l*b)/p (m*g*l*(M+m))/p     0];
         
B=[     0;                                % input matrix
    (I+m*l^2)/p;      
        0; 
    (m*l)/p];

C=[1 0 0 0;                               % output matrix
   0 0 0 1];

D=zeros(2,1);                             % feedforward matrix

states={'x' 'xdot' 'phi' 'phidot'};
inputs={'u'};
outputs={'x' 'phidot'};

sys_ss=ss(A,B,C,D,'statename', states, 'inputname', inputs, 'outputname', outputs);

%% Linear Quadratic Regulator Design

ct=rank(ctrb(A,B));  % test for controllability

Q=[50 0 0 0;         % penalty on x
   0 0 0 0;
   0 0 25 0;         % penalty on phi
   0 0 0 0];  

R=[1];               % penalty on u

K=lqr(A,B,Q,R);      % LQR gain matrix

Ac=(A-(B*K));        % modified system matrix


%% Precompensation to reduce steady state error

N=-1/(C*inv(Ac)*B);  % scaling factor calculation
Sf=N(1,1);             

sys_cl=ss(Ac,Sf*B,C,D,'statename', states, 'inputname', inputs, 'outputname', outputs);

%% Estimator Design

ot=rank(obsv(A,C));             % test for observability
P=[-40 -41 -42 -43];            % desired estimator poles
L=place(A',C',P)';              % estimator gain

Ace=[     Ac          B*K;      % modified system matrix
     zeros(size(A)) (A-L*C)];   
Bce=[Sf*B; zeros(size(B))];     % modified input matrix
Cce=[C zeros(size(C))];         % modified output matrix
Dce=0;                          % feedforward matrix

%% Kalman Filter

Qk=0.005;          % process noise covariance (motion model)
Rk=0.001;          % measurement noise covariance (sensor model)
