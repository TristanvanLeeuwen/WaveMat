classdef opGrad < opSpot
%OpGrad n-dimensional Gradient operaror based on spectral
%differences. Note that Div = -Grad^T !
%
% use:
%   Grad = opGrad(N,L,method,store)
%
% input:
%   N      - number of gridoints in each direction
%   L      - size of physical domain in m
%   method - 'fourier' or 'cheb'
%   store  - whether to store or generate on-the-fly

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Properties
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   properties ( Access = private )
      funHandle       % Multiplication function
   end % Properties
   
   properties ( SetAccess = private, GetAccess = public )
      A               % array for matrices
      x
      N
      L
      method
      store
   end % properties
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - public
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % opGrad Constructor.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function op = opGrad(N,L,method,store)
         % number of dimensions
         nd = length(N);
         
         op = op@opSpot('Grad',nd*prod(N),prod(N));
         op.A           = cell(1,nd);
         op.x           = cell(1,nd);
         op.N           = N;
         op.L           = L;
         op.method      = method;
         op.store       = store;
         op.cflag       = false;
         op.sweepflag   = false;
         
         % initialize
         switch method
             case 'fourier'
                 for i = 1:nd
                     op.x{i} = 2*pi*(0:N(i)-1)'/N(i);    
                     a       = L(i)/(op.x{i}(end) - op.x{i}(1));
                     op.x{i} = a*(op.x{i} - op.x{i}(1));
                     if store
                        [~,op.A{i}] = fourdif(N(i),1);
                        op.A{i}     = (1/a)*op.A{i};
                     else
                        op.A{i} = @(x,mode)((1/a)*fourdifft(x,1,mode));
                     end
                 end
             case 'cheb'
                 for i = 1:nd
                     op.x{i} = sin(pi*[N(i)-1:-2:1-N(i)]'/(2*(N(i)-1)));
                     a       = L(i)/(op.x{i}(end) - op.x{i}(1));
                     op.x{i} = a*(op.x{i} - op.x{i}(1));
                     if store
                         [~,op.A{i}] = chebdif(N(i),1);
                         op.A{i}     = (1/a)*op.A{i};
                     else
                         op.A{i} = @(x,mode)((1/a)*chebdifft(x,1,mode));
                     end
                 end
             case 'fd'
                 for i = 1:nd
                     op.x{i} = linspace(0,1,N(i));
                     a       = L(i)/(op.x{i}(end) - op.x{i}(1));
                     op.x{i} = a*(op.x{i} - op.x{i}(1));
                     if store
                         [~,op.A{i}] = fddif(N(i),1);
                         op.A{i}     = (1/a)*op.A{i};
                     else
                         op.A{i} = @(x,mode)((1/a)*fddif(N(i),1)*x);
                     end
                 end
                 
             otherwise
                 error('Method not supported...');
         end
         
         
         % Create function handle
         op.funHandle = @opGrad_multiply;
      end % function 
      
   end % methods - public
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - protected
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Access = protected )
      
      % Multiplication
      function y = multiply(op,x,mode)
         y = op.funHandle(op,x,mode);
      end % Multiply
      
   end % Methods
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - private
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Access = private )
      
       function y = opGrad_multiply(op,x,mode)
           nd    = length(op.N);
           if mode == 1
               y = zeros(prod(op.N),nd);
               for j=1:nd
                   t = x;
                   for i=1:nd
                       k = op.N(i);
                       t = reshape(t,k,[]);
                       if i==j
                           if op.store
                               z = op.A{i}' * t;
                           else
                               z = op.A{i}(t,-1);
                           end
                       else
                           z = t;
                       end
                       t = z.';
                   end
                   y(:,j) = t(:);
               end
               y = y(:);
           else
               x = reshape(x,[prod(op.N) nd]);
               y = zeros(prod(op.N),1);
               for j=1:nd
                   t = x(:,j);
                   for i=1:nd
                       k = op.N(i);
                       t = reshape(t,k,[]);
                       if i == j
                           if op.store
                               z = op.A{i} * t;
                           else
                               z = op.A{i}(t,1);
                           end
                       else
                           z = t;
                       end
                       t = z.';
                   end
                   y = y + t(:);
               end
           end
           
       end
      
   end % Methods
   
end % classdef

