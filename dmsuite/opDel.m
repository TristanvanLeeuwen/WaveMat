classdef opDel < opSpot
%OpDel 1st-order derivative operator based on spectral differences
%
% use:
%   Del = opDel(N,L,dir,method,store)
%
% input:
%   N      - number of gridoints in each direction
%   L      - size of physical domain in m
%   dir    - index indicating which direction to act on
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
      dir
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
      function op = opDel(N,L,dir,method,store)
         % number of dimensions
         nd = length(N);
         
         op = op@opSpot('Del',prod(N),prod(N));
         op.A           = cell(1,nd);
         op.x           = cell(1,nd);
         op.N           = N;
         op.L           = L;
         op.dir         = dir;
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
                     if i==op.dir
                         if store
                            [~,op.A] = fourdif(N(i),1);
                            op.A     = (1/a)*op.A;
                         else
                            op.A = @(x,mode)((1/a)*fourdifft(x,1,mode));
                         end
                     end
                 end
             case 'cheb'
                 for i = 1:nd
                     op.x{i} = sin(pi*[N(i)-1:-2:1-N(i)]'/(2*(N(i)-1)));
                     a       = L(i)/(op.x{i}(end) - op.x{i}(1));
                     op.x{i} = a*(op.x{i} - op.x{i}(1));
                     if i==dir
                         if store
                             [~,op.A] = chebdif(N(i),1);
                             op.A     = (1/a)*op.A;
                         else
                             op.A = @(x,mode)((1/a)*chebdifft(x,1,mode));
                         end
                     end
                 end
             otherwise
                 error('Method not supported...');
         end
         
         
         % Create function handle
         op.funHandle = @opDel_multiply;
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
      
       function y = opDel_multiply(op,x,mode)
           nd    = length(op.N);
           if mode == 1
               for i=1:nd
                   x = reshape(x,op.N(i),[]);
                   if i==op.dir
                       if op.store
                           z = op.A * x;
                       else
                           z = op.A(x,1);
                       end
                   else
                       z = x;
                   end
                   x = z.';
               end
               y = x(:);
           else
               for i=1:nd
                   x = reshape(x,op.N(i),[]);
                   if i==op.dir
                       if op.store
                           z = op.A' * x;
                       else
                           z = op.A(x,-1);
                       end
                   else
                       z = x;
                   end
                   x = z.';
               end
               y = x(:);
           end
           
       end
      
   end % Methods
   
end % classdef

