clc
clear all
close all
[userview systemview] = memory
%%
% Pride equations ( <http//geodus1.ta.tudelft.nl/PrivatePages/C.P.A.Wapenaar/4_Journals/J.Appl.Mech/AppM_04.pdf> )
%
% $$\rho_b \frac{\partial \bf{v}_s}{\partial t}+\rho_f\frac{\partial \bf{w}}{\partial t}$$
%
% $$\epsilon \frac{\partial \bf{E}}{\partial t} - \nabla \times \bf{H} +\sigma \bf{E}=-\bf{J}^e$$
%
% $$\mu \frac{\partial \bf{H}}{\partial t} + \nabla \times \bf{E}= -\bf{J}^m$$
%
% Reciprocity theorems for diffusion, flow and waves, Wapenaar, Fokkema
% Matrix form
%
% $$j\omega \left[\begin{array}{ccccccccccccc} \rho_b & \bf{0} & \bf{0} & \bf{0}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \rho_f & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0} & \rho_b & \bf{0} & \bf{0}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \rho_f & \bf{0} & \bf{0} \\
%   \bf{0} & \bf{0} & \rho_b & \bf{0}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \rho_f & \bf{0} \\
%   \bf{0} & \bf{0} & \bf{0} & \bf{I}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0} & \bf{0} & \bf{0} & \bf{0}    & \bf{I}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0} & \bf{0} & \bf{0} & \bf{0}    & \bf{0}  & \bf{I}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0} & \bf{0} & \bf{0} & \bf{0}    & \bf{0}  & \bf{0}  & \bf{I} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0} & \bf{0} & \bf{0} & \bf{0}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{I} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0} & \bf{0} & \bf{0} & \bf{0}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{I} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \rho_f & \bf{0} & \bf{0} & \bf{0}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0} & \rho_f & \bf{0} & \bf{0}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0} & \bf{0} & \rho_f & \bf{0}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0} & \bf{0} & \bf{0} & \bf{0}    & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{I} \\
%   \end{array}\right]\bf{u}
%   \left[\begin{array}{ccccccccccccc} 
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \frac{\eta}{k} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \frac{\eta}{k} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \frac{\eta}{k} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \end{array}\right]\bf{u}+
%   \left[\begin{array}{ccccccccccccc} 
%   \bf{I}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{I} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{0} & \bf{I} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{0} & \bf{0} & {\lambda+2\mu} & \lambda & \lambda & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{d1} \\
%   \bf{0}& \bf{0} & \bf{0} & \lambda & {\lambda+2\mu} & \lambda & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{d1} \\
%   \bf{0}& \bf{0} & \bf{0} & \lambda & \lambda & {\lambda+2\mu} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{d1} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \mu    & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \mu    & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \mu    & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{I} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{I} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{I} & \bf{0} \\ 
%   \bf{0}& \bf{0} & \bf{0} & \bf{d1}        & \bf{d1} & \bf{d1} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & M\bf{I}\\ 
%   \end{array}\right]\left[\begin{array}{ccccccccccccc}
%   \bf{0}& \bf{0} & \bf{0} & \bf{L1}        & \bf{0}  & \bf{0}  & \bf{0} & \bf{L3 & \bf{L2}& \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{L2} & \bf{0}  & \bf{L3}& \bf{0} & \bf{L1}& \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{L3} & \bf{L2}& \bf{L1}& \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{L1}&\bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{L2}& \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{0} & \bf{L3}& \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{0}& \bf{L3}& \bf{L2}& \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{L3}&\bf{0} & \bf{L1}& \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\ 
%   \bf{L2}&\bf{L1}& \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} \\
%   \bf{0} &\bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{L1}\\
%   \bf{0} &\bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{L2}\\
%   \bf{0} &\bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{0} & \bf{L3}\\
%   \bf{0} &\bf{0} & \bf{0} & \bf{0}         & \bf{0}  & \bf{0}  & \bf{0} & \bf{0} & \bf{0} & \bf{L1}& \bf{L2}& \bf{L3}& \bf{0} \\
% \end{array}\right]\bf{u}=\bf{s}$$
%
% where
%
% $$\left[\begin{array}{c} -\bf{J}^e\\ -\bf{J}^m\\ \end{array}\right]$$
%
% and
%
% $$\left[\begin{array}{ccc} \bf{0} & -\frac{\partial}{\partial z} & \frac{\partial}{\partial y} \\  \frac{\partial}{\partial z} & \bf{0} & -\frac{\partial}{\partial x} \\ -\frac{\partial}{\partial y}& \frac{\partial}{\partial x}& \bf{0}\\
% \end{array}\right]$$

N        = 20            ;
[x,D]    = fourdif(N,1)   ;  
I        = speye(N)       ;
L1       = kron(I,D)      ;
L2       = kron(D,I)      ;
L3       = sparse(N^2,N^2) ;
L0       = sparse(N^2,N^2) ;
y        = x              ;
tspan    = 0:0.1:30      ;
[X,Y]    = meshgrid(x)    ;
X        = (X)      ;
Y        = (Y)      ;

eta      = .9;
k        = speye(3*N^2);
L        = speye(3*N^2);%zeros(3*N^2,3*N^2);
rho_B    = 2.6e3;
rho_b    = rho_B*speye(3*N^2);
rho_F    = 1e3;
rho_f    = rho_F*speye(3*N^2);
Kfr      = .1;
Kf       = 1;
Ks       = 1;
phi      =.1;
NN       = 1;
delta    = Kf/(phi*Ks^2)*((1-phi)*Ks-Kfr);
Kg       = (Kfr+phi*Kf+(1+phi)*Ks*delta)/(1+delta);
dd       = (Kf+Ks*delta)/(1+delta);
MM       = (1/phi*(Kf/(1+delta)));
VV       = (Kg-2/3*NN);

% D matrix
L5122a = sparse([L1, zeros(N^2,N^2), zeros(N^2,N^2); 
                 zeros(N^2,N^2), L2, zeros(N^2,N^2); 
                 zeros(N^2,N^2), zeros(N^2,N^2), L3]);
L5122b = sparse([...
           zeros(N^2,N^2)    ,  L3   ,  L2;
           (L3)  ,  zeros(N^2,N^2)   ,  L1;
           (L2)  ,  (L1) ,  zeros(N^2,N^2)]);

L5122 = sparse([...
            zeros(3*N^2,3*N^2)  , L5122a            , L5122b;
            -transpose(L5122a)  , zeros(3*N^2,3*N^2), zeros(3*N^2,3*N^2);
            -transpose(L5122b)  , zeros(3*N^2,3*N^2), zeros(3*N^2,3*N^2)]);

L5133 = sparse([zeros(N^2,N^2), zeros(N^2,N^2) , zeros(N^2,N^2) , L1;
                zeros(N^2,N^2), zeros(N^2,N^2) , zeros(N^2,N^2) , L2;
                zeros(N^2,N^2), zeros(N^2,N^2) , zeros(N^2,N^2) , L3;
                -transpose(L1), -transpose(L2) , -transpose(L3) ,zeros(N^2,N^2)]);

L51   = sparse([...
             L5122 ,  zeros(9*N^2,4*N^2);
             zeros(4*N^2,9*N^2) , L5133]);

% C matrix
% L222a = sparse([...      %3D
%                 [VV+2*NN]*eye(N^2), [VV*NN]*eye(N^2) , [VV*NN]*eye(N^2) ;
%                 (VV*NN*eye(N^2))  , [VV+2*NN]*eye(N^2), [VV*NN]*eye(N^2) ;
%                 (VV*NN*eye(N^2)  ), (VV*NN*eye(N^2))  , [VV+2*NN]*eye(N^2) ]);
L222a = sparse([...   %2D
                [VV+2*NN]*eye(N^2)     , [VV]*eye(N^2) , zeros(N^2,N^2) ;
                (VV*eye(N^2)) , [VV+2*NN]*eye(N^2)     , zeros(N^2,N^2) ;
                zeros(N^2,N^2), zeros(N^2,N^2)         , eye(N^2,N^2) ]);            
L222b = sparse([...   %2D
                zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
                zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
                zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2)]);
% L222c = sparse([...   %3D
%                 [NN]*eye(N^2), zeros(N^2,N^2),zeros(N^2,N^2);%
%                 zeros(N^2,N^2), [NN]*eye(N^2),zeros(N^2,N^2);%
%                 zeros(N^2,N^2), zeros(N^2,N^2),[NN]*eye(N^2)]);
L222c = sparse([...
                eye(N^2), zeros(N^2,N^2),zeros(N^2,N^2);%
                zeros(N^2,N^2), eye(N^2),zeros(N^2,N^2);%
                zeros(N^2,N^2), zeros(N^2,N^2),[NN]*eye(N^2)]);

            % 
L222 = sparse([...
                eye(3*N^2) ,zeros(3*N^2,3*N^2) ,zeros(3*N^2,3*N^2);
                zeros(3*N^2,3*N^2), L222a , L222b; %
                zeros(3*N^2,3*N^2),(L222b), L222c]);%
L223  = sparse([...
             zeros(1*N^2,3*N^2),zeros(1*N^2,1*N^2);
             zeros(1*N^2,3*N^2),zeros(1*N^2,1*N^2);
             zeros(1*N^2,3*N^2),zeros(1*N^2,1*N^2);
             zeros(1*N^2,3*N^2),dd*eye(N^2);%zeros(1*N^2,1*N^2);%
             zeros(1*N^2,3*N^2),dd*eye(N^2);%zeros(1*N^2,1*N^2);%
             zeros(1*N^2,3*N^2),zeros(1*N^2,1*N^2);%%%%dd*eye(N^2);%
             zeros(1*N^2,3*N^2),zeros(1*N^2,1*N^2);%%%%dd*eye(N^2);%
             zeros(1*N^2,3*N^2),zeros(1*N^2,1*N^2);%%%%dd*eye(N^2);%
             zeros(1*N^2,3*N^2),zeros(1*N^2,1*N^2)]);%dd*eye(N^2)]);%
L233  = sparse([...
                eye(3*N^2,3*N^2) , zeros(3*N^2,1*N^2);
                zeros(1*N^2,3*N^2), MM*eye(1*N^2,1*N^2)]);
L22 =  sparse([...
                L222  ,L223;%
                transpose(L223) ,L233]);%
                
% B matrix  
L5233  = sparse([eta.*inv(k)   ,   zeros(3*N^2,N^2);    zeros(N^2,3*N^2), zeros(N^2,N^2)]);
L5222  = sparse(zeros(9*N^2,9*N^2)) ;


L522213 = sparse(zeros(9*N^2,4*N^2));
             
L52   = sparse([...
                L5222 , L522213; 
                transpose(L522213) , L5233 ]);
% A matrix             
% M22      = sparse([rho_b             ,zeros(3*N^2,3*N^2),zeros(3*N^2,3*N^2);
%                    zeros(3*N^2,3*N^2),eye(3*N^2)        ,zeros(3*N^2,3*N^2);
%                    zeros(3*N^2,3*N^2),zeros(3*N^2,3*N^2),eye(3*N^2)]);     
% M23      = sparse([rho_f             ,zeros(3*N^2,1*N^2);
%                    zeros(3*N^2,3*N^2),zeros(3*N^2,1*N^2);
%                    zeros(3*N^2,3*N^2),zeros(3*N^2,1*N^2)]);
% M33      = sparse([zeros(3*N^2,3*N^2),zeros(3*N^2,1*N^2);
%                    zeros(1*N^2,3*N^2),eye(N^2)]);
% % M22      = sparse([rho_B*eye(N^2)    ,zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),eye(N^2)    ,zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),eye(N^2)          ,zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),eye(1*N^2)        ,zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),eye(1*N^2)        ,zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),eye(1*N^2)        ,zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),eye(1*N^2)        ,zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),eye(1*N^2)        ,zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),eye(1*N^2)        ]);
% % M23      = sparse([rho_F*eye(N^2)    ,zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2)]);
% % M33      = sparse([zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2);
% %                    zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),zeros(1*N^2,1*N^2),eye(N^2)         ]);
% % 
% % M        = sparse([  M22           , M23;
% %                     transpose(M23) , M33]);


%
%QQQ = abs(M\speye(13*N^2));
M  = sparse([...
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),1e-3*eye(N^2) ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),1e-3*eye(N^2) ,zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),eye(N^2)      ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),eye(N^2)      ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),eye(N^2)      ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),eye(N^2)      ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),eye(N^2)      ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),eye(N^2)      ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),eye(N^2)      ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
             1e-3*eye(N^2) ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),-eye(N^2)*1e-4,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),1e-3*eye(N^2) ,zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),-eye(N^2)*1e-4,zeros(N^2,N^2),zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),eye(N^2)      ,zeros(N^2,N^2);
             zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),zeros(N^2,N^2),eye(N^2)]);
L5 = sparse(L22*L51+L52); 

%init cond
Vsx    = sparse(N,N);
Vsy    = sparse(N,N);
Vsz    = sparse(N,N);
Tau_xx = sparse(N,N);
Tau_yy = sparse(N,N);
Tau_zz = sparse(N,N);
Tau_yz = sparse(N,N);
Tau_zx = sparse(N,N);
Tau_xy = sparse(N,N);
Wx     = sparse(N,N);
Wy     = sparse(N,N);
Wz     = sparse(N,N);
P      = sparse(exp(-10*((X - mean(x)).^2 + (Y - mean(x)).^2)));


vsx    = reshape(Vsx,N^2,1)     ;
vsy    = reshape(Vsy,N^2,1)     ;
vsz    = reshape(Vsz,N^2,1)     ;
tau_xx = reshape(Tau_xx,N^2,1)  ;
tau_yy = reshape(Tau_yy,N^2,1)  ;
tau_zz = reshape(Tau_zz,N^2,1)  ;
tau_yz = reshape(Tau_yz,N^2,1)  ;
tau_zx = reshape(Tau_zx,N^2,1)  ;
tau_xy = reshape(Tau_xy,N^2,1)  ;
wx     = reshape(Wx,N^2,1)      ;
wy     = reshape(Wy,N^2,1)      ;
wz     = reshape(Wz,N^2,1)      ;
p      = reshape(P,N^2,1)      ;

options  = odeset('OutputFcn',@odewbar),
[t,usol] = ode113(@(t,uu)-(L5*uu),[0 max(tspan)],[vsx;vsy;vsz;tau_xx;tau_yy;tau_zz;tau_yz;tau_zx;tau_xy;wx;wy;wz;p],options) ;
% vidObj = VideoWriter('B=0.avi');
% vidObj.FrameRate=length(t);
% open(vidObj);
for j=1:length(t);
    uplot1 =reshape(usol(j,1:N^2),N,N);
    uplot2 =reshape(usol(j,N^2+1:2*N^2),N,N);
    uplot3 =reshape(usol(j,2*N^2+1:3*N^2),N,N);
    uplot4 =reshape(usol(j,3*N^2+1:4*N^2),N,N);
    uplot5 =reshape(usol(j,4*N^2+1:5*N^2),N,N);
    uplot6 =reshape(usol(j,5*N^2+1:6*N^2),N,N);
    uplot7 =reshape(usol(j,6*N^2+1:7*N^2),N,N);
    uplot8 =reshape(usol(j,7*N^2+1:8*N^2),N,N);
    uplot9 =reshape(usol(j,8*N^2+1:9*N^2),N,N);
    uplot10=reshape(usol(j,9*N^2+1:10*N^2),N,N);
    uplot11=reshape(usol(j,10*N^2+1:11*N^2),N,N);
    uplot12=reshape(usol(j,11*N^2+1:12*N^2),N,N);
    uplot13=reshape(usol(j,12*N^2+1:13*N^2),N,N);

subplot(5,4,1)
    imagesc(x,x,uplot1,[-1 1]*.1);title('V^s_x');colormap(gray);axis equal tight;
subplot(5,4,2)
    imagesc(x,x,uplot2,[-1 1]*.1);title('V^s_y');colormap(gray);axis equal tight;
subplot(5,4,3)
    imagesc(x,x,uplot3,[-1 1]*.1);title('V^s_z');colormap(gray);axis equal tight;
subplot(5,4,4)
    imagesc(x,x,uplot4,[-1 1]*.1);title('\tau_{xx}');colormap(gray);axis equal tight;
subplot(5,4,5)
    imagesc(x,x,uplot5,[-1 1]*.1);title('\tau_{yy}');colormap(gray);axis equal tight;
subplot(5,4,6)
    imagesc(x,x,uplot6,[-1 1]*.1);title('\tau_{zz}');colormap(gray);axis equal tight;
subplot(5,4,7)
    imagesc(x,x,uplot7,[-1 1]*.1);title('\tau_{yz}');colormap(gray);axis equal tight;
subplot(5,4,8)
    imagesc(x,x,uplot8,[-1 1]*.1);title('\tau_{zx}');colormap(gray);axis equal tight;
subplot(5,4,9)
    imagesc(x,x,uplot9,[-1 1]*.1);title('\tau_{xy}');colormap(gray);axis equal tight;
subplot(5,4,10)
    imagesc(x,x,uplot10,[-1 1]*.1);title('w_x');colormap(gray);axis equal tight;
subplot(5,4,11)
    imagesc(x,x,uplot11,[-1 1]*.1);title('w_y');colormap(gray);axis equal tight;
subplot(5,4,12)
    imagesc(x,x,uplot12,[-1 1]*.1);title('w_z');colormap(gray);axis equal tight; 
subplot(5,4,13)
    imagesc(x,x,uplot13,[-1 1]*.1);title('p');colormap(gray);axis equal tight;    
    drawnow
    pause(0.01)

%       F=getframe(gcf);
%       writeVideo(vidObj,F);

end
% close(vidObj);
%%
% 
%  
% Reconsider Maxwell equations, by
% multiplying the first and second equation with $\nabla$,
% and $\frac{\partial}{\partial t}$ respectively
%
% This results in one combined equation (the electric field equation)
%
% $$\epsilon \frac{\partial^2 \bf{E}}{\partial t^2}+\sigma
% \frac{\partial \bf{E}}{\partial t}+\nabla \times 
% (\frac{1}{\mu} \nabla \times \bf{E})=-\frac{\partial \bf{J}^e}{\partial t}-
% \nabla \times (\frac{1}{\mu} \bf{J}^m)$$
%
% This results in one combined equation (the magnetic field equation)
%
% $$\mu \frac{\partial^2 \bf{H}}{\partial t^2}+\mu
% \frac{\partial \bf{H}}{\partial t}+\nabla \times 
% (\frac{1}{\sigma} \nabla \times \bf{H})=-\bf{J}^m +
% \nabla \times (\frac{1}{\sigma} \bf{J}^e)$$


%%
% Appendix, first put fourdif functionin seperate file to run file above
%%

%     function [x, DM] = fourdif(N,m)
% %
% % The function [x, DM] = fourdif(N,m) computes the m'th derivative Fourier 
% % spectral differentiation matrix on grid with N equispaced points in [0,2pi)
% % 
% %  Input:
% %  N:        Size of differentiation matrix.
% %  M:        Derivative required (non-negative integer)
% %
% %  Output:
% %  x:        Equispaced points 0, 2pi/N, 4pi/N, ... , (N-1)2pi/N
% %  DM:       m'th order differentiation matrix
% %
% % 
% %  Explicit formulas are used to compute the matrices for m=1 and 2. 
% %  A discrete Fouier approach is employed for m>2. The program 
% %  computes the first column and first row and then uses the 
% %  toeplitz command to create the matrix.
% 
% %  For m=1 and 2 the code implements a "flipping trick" to
% %  improve accuracy suggested by W. Don and A. Solomonoff in 
% %  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
% %  The flipping trick is necesary since sin t can be computed to high
% %  relative precision when t is small whereas sin (pi-t) cannot.
% %
% %  S.C. Reddy, J.A.C. Weideman 1998.  Corrected for MATLAB R13 
% %  by JACW, April 2003.
%  
% 
%     x=2*pi*(0:N-1)'/N;                       % gridpoints
%     h=2*pi/N;                                % grid spacing
%     zi=sqrt(-1);
%     kk=(1:N-1)';
%     n1=floor((N-1)/2); n2=ceil((N-1)/2);
%     if m==0,                                 % compute first column
%       col1=[1; zeros(N-1,1)];                % of zeroth derivative
%       row1=col1;                             % matrix, which is identity
% 
%     elseif m==1,                             % compute first column
%       if rem(N,2)==0                         % of 1st derivative matrix
% 	topc=cot((1:n2)'*h/2);
%         col1=[0; 0.5*((-1).^kk).*[topc; -flipud(topc(1:n1))]]; 
%       else
% 	topc=csc((1:n2)'*h/2);
%         col1=[0; 0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];
%       end;
%       row1=-col1;                            % first row
% 
%     elseif m==2,                             % compute first column  
%       if rem(N,2)==0                         % of 2nd derivative matrix
% 	topc=csc((1:n2)'*h/2).^2;
%         col1=[-pi^2/3/h^2-1/6; -0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];
%       else
% 	topc=csc((1:n2)'*h/2).*cot((1:n2)'*h/2);
%         col1=[-pi^2/3/h^2+1/12; -0.5*((-1).^kk).*[topc; -flipud(topc(1:n1))]];
%       end;
%       row1=col1;                             % first row 
% 
%     else                                     % employ FFT to compute
%       N1=floor((N-1)/2);                     % 1st column of matrix for m>2
%       N2 = (-N/2)*rem(m+1,2)*ones(rem(N+1,2));  
%       mwave=zi*[(0:N1) N2 (-N1:-1)];
%       col1=real(ifft((mwave.^m).*fft([1 zeros(1,N-1)])));
%       if rem(m,2)==0,
% 	row1=col1;                           % first row even derivative
%       else
% 	col1=[0 col1(2:N)]'; 
% 	row1=-col1;                          % first row odd derivative
%       end;
%     end;
%     DM=toeplitz(col1,row1);                   


%%
% 
%  
% 
%  
%% References:
% R.J.Leveque, Finite Volume methods for hyperbolic problems, http://csclub.uwaterloo.ca/~lbovard/finite-volume/Finite_Volume_Methods_for_Hyperbolic_Problems.pdf
%
% S.Phadke, S.Yerneni, A PVM Implementation of 2D Acoustic Wave Modelling on PARAM 10000, http://www.cdac.in/html/pdf/pvmimpl1.pdf
%
% L.N.Tretethen, Spectral Methods in Matlab,
% http://www.mathworks.ir/downloads/Spectral%20Methods%20in%20MATLAB%5Bwww.mathworks.ir%5D.pdf 
% http://people.maths.ox.ac.uk/trefethen/spectral.html
%
% J.N.Kutz, AMATH 581 Practical Scientific Computing, http://courses.washington.edu/amath581/581.pdf
% http://depts.washington.edu/amath/courses/571-winter-2003/matlab.html
% http://github.com/pycckuu/ScientificComputation/tree/c9fff18f8ebb4cd8d8ecdd98450b7127b201485b
% 
% G.Chen, B.Cloutier, N.Li, B.K.Muite and P.Rigge, Parallel Spectral Numerical Methods,
% http://shodor.org/media/content//petascale/materials/UPModules/Parallel_Spectral_Methods/ParallelNumericalMethods.pdf
%
% R.J. Leveque, Finite Volume Methods, http://csclub.uwaterloo.ca/~lbovard/finite-volume/Finite_Volume_Methods_for_Hyperbolic_Problems.pdf
%
% GNU Octave, Memory efficiency, http://www.gnu.org/software/octave
% http://www.gnu.org/software/octave/doc/interpreter/Creating-Sparse-Matrices.html
%
% Mathworks, Memory optimization, http://www.mathworks.nl/help/matlab/numeric-types.html
% http://blogs.mathworks.com/loren/2012/02/06/using-gpus-in-matlab/
% http://www.mathworks.nl/help/matlab/matlab_prog/profiling-for-improving-performance.html
%
% http://www.cis.upenn.edu/~cis515/cis515-11-sl9.pdf\
% http://www.geophysicist.nl/lib/exe/fetch.php?media=wapenaareage2007abs.pdf
%
% http://stackoverflow.com/questions/18278334/how-to-make-a-movie-of-a-subplot-in-which-each-subplot-is-a-frame
% 
M=   [ 1e3 0 0 0 0 0 0 1e4 0 0;
       0 1e3 0 0 0 0 0 0 1e4 0;
       0 0 1 0 0 0 0 0 0 0;
       0 0 0 1 0 0 0 0 0 0;
       0 0 0 0 1 0 0 0 0 0;
       0 0 0 0 0 1 0 0 0 0;
       0 0 0 0 0 0 1 0 0 0;
       1e4 0 0 0 0 0 0 0 0 0;
       0 1e4 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 1];

M=   [ 1e3 0 0 0 0 0 0 1e4 0 0;
       0 1 0 0 0 0 0 0 0 0;
       0 0 1 0 0 0 0 0 0 0;
       0 0 0 1 0 0 0 0 0 0;
       0 0 0 0 1 0 0 0 0 0;
       0 0 0 0 0 1 0 0 0 0;
       0 0 0 0 0 0 1 0 0 0;
       1e4 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 1 0;
       0 0 0 0 0 0 0 0 0 1];   