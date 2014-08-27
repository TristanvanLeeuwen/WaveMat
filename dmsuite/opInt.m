classdef opInt < opSpot
% opInt - Interpolation based on Chebychev or Fourier
%


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Properties
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   properties ( Access = private )
      funHandle       % Multiplication function
      xc
      x
      method
   end % Properties
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - public
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % opDFT. Constructor.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function op = opInt(xc,x,method)
         
         op = op@opSpot('Int',length(x),length(xc));
         op.cflag       = false;
         op.sweepflag   = true;
         op.x           = x(:);
         op.xc          = xc(:);
         
         % Create function handle
         switch method
             case 'cheb'
                 op.funHandle = @opInt_cheb;
             case 'fourier'
                 op.funHandle = @opInt_fourier;
             case 'fd'
                 op.funHandle = @opInt_fd;  
             otherwise
                 error('Method not supported');
         end

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
      
      function y = opInt_cheb(op,f,mode) 
          % map grid to [-1,1]
          a = 2./(op.xc(end) - op.xc(1));
          b = -1 -op.xc(1).*a;
          % interpolate
          y = chebint(f,a*op.x + b,mode); 
      end
      
      function y = opInt_fourier(op,f,mode)  
          % map grid to [0,2pi]
          a = 2*pi./(op.xc(end) - op.xc(1));
          b =-op.xc(1).*a;
          
          y = fourint(f,a*op.x + b,mode); 
      end
      
      function y = opInt_fd(op,f,mode)
          y = f;
      end
      
   end % Methods
   
end % classdef

