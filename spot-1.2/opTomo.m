% -----------------------------------------------------------------------
% Copyright 2013 Centrum Wiskunde & Informatica, Amsterdam
%
% Author:  Folkert Bleichrodt
% Contact: F.Bleichrodt@cwi.nl
% 
% This file is the ASTRA-Spot operator "opTomo" to be used with the
% All Scale Tomographic Reconstruction Antwerp Toolbox ("ASTRA Toolbox").
%
% The ASTRA-Spot operator is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The ASTRA-Spot operator is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the ASTRA-Spot operator. If not, see <http://www.gnu.org/licenses/>.
%
%-----------------------------------------------------------------------

classdef opTomo < opSpot
    % opTomo - Wrapper for ASTRA tomography projector


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = private )
        % Multiplication function
        funHandle
        proj_geom
        vol_geom
        proj_handle
    end % Properties

    properties ( SetAccess = private, GetAccess = public )
        proj_size
        vol_size
    end % properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % opTomo. Constructor.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opTomo(type, proj_geom, vol_geom)
            % determine the dimension
            % TODO: Check consistency between proj_geom and vol_geom
            is2D = ~isfield(vol_geom, 'GridSliceCount');

            if is2D
                proj_id = astra_create_projector(type, proj_geom, vol_geom);
                proj    = astra_projector_handle(proj_id);
            else
                if ~strcmp(type, 'cuda')
                    error(['Only type ' 39 'cuda' 39 ' is supported ' ...
                           'for 3D geometries.'])
                end
                proj = [];
                if strcmp(proj_geom.type, 'parallel3d_vec') || ...
                    strcmp(proj_geom.type, 'cone_vec')
                    angleCount = size(proj_geom.Vectors, 1);
                else
                    angleCount = numel(proj_geom.ProjectionAngles);
                end
            end

            proj_size = astra_geom_size(proj_geom);
            vol_size  = astra_geom_size(vol_geom);

            % construct operator
            op = op@opSpot('opTomo', prod(proj_size), prod(vol_size));

            op.proj_size = proj_size;
            op.vol_size  = vol_size;
            
            % Wrap the projector in a handle class to allow a destructor
            op.proj_handle = proj;
            op.proj_geom   = proj_geom;
            op.vol_geom    = vol_geom;

            op.cflag       = false;
            op.sweepflag   = false;

            % Create function handle
            if is2D
                op.funHandle = @opTomo_intrnl2D;
            else
                op.funHandle = @opTomo_intrnl3D;
            end
        end %

    end % methods - public

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )

        % Multiplication
        function y = multiply(op,x,mode)
            y = op.funHandle(op, x, mode);
        end % Multiply

    end % Methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - private
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = private )

        % 2D projection code
        function y = opTomo_intrnl2D(op,x,mode)
            if mode == 1
                % x is passed as a vector. Reshape it into an image.
                x = reshape(x, op.vol_size);

                % forward projection
                [sid, s] = astra_create_sino(x, op.proj_handle.id);
                astra_mex_data2d('delete', sid);

                % Now s is the sinogram. Reshape it back into a vector
                y = s(:);
            else
                % x is passed as a vector. Reshape it into a sinogram.
                x = reshape(x, op.proj_size);

                % backprojection
                [vid,v] = astra_create_backprojection(x, op.proj_handle.id);
                astra_mex_data2d('delete', vid);

                % Now v is the resulting volume. Reshape it back into a vector
                y = v(:);
            end
        end


        % 3D projection code
        function y = opTomo_intrnl3D(op,x,mode)
            if mode == 1
                % X is passed as a vector. Reshape it into a volume.
                x = reshape(x, op.vol_size);

                % forward projection
                [sid, s] = astra_create_sino3d_cuda(x, op.proj_geom, op.vol_geom);
                astra_mex_data3d('delete', sid);

                % now s is the projection data. Reshape it back into a vector
                y = s(:);
            else
                % X is passed as a vector. Reshape it into projection data.
                x = reshape(x, op.proj_size);

                % forward projection
                [vid, v] = astra_create_backprojection3d_cuda(x, op.proj_geom, op.vol_geom);
                astra_mex_data3d('delete', vid);

                % Now v is the resulting volume. Reshape it back into a vector
                y = v(:);
            end
        end

    end % Methods

end % classdef
