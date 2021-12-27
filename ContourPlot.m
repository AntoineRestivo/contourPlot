classdef ContourPlot < handle
% Describes an on-going contour plot
%
% The plot area is divided into a grid of ``div(1)+1`` x ``div(2)+1`` points. For the x axis, there is a linear mapping such that
%
% ``0 -> xmin``
% ``div(1)+1 -> xmax``
%
% and the same for the y axis:
%
% ``0 -> ymin``
% ``div(2)+1 -> ymax``
%
% The plot area is encoded into `.axis`.
%
% We can take ``div(1)``, ``div(2)`` a pretty large power-of-two, knowing that the base memory usage is proportional to ``div(2)``. Thus, something
% like ``2^16`` is optimal. The ``div`` corresponds to the finest refinement possible, when our code supports
% refinement.
%
% We start the marching squares algorithm with a step size of ``s``. By default, this is ``2^8``.
%
% The method `.evaluate` takes two integer coordinates, and either computes the function value if it has not
% been computed, or returns the value from the cache. Values are cached in a sparse matrix `.values`
% shifting the indices by 1 (because our integer indices are 0-based, and Matlab is 1-based). In the matrix, the
% value ``0`` is not stored as is, but using a value ``eps``, which is approximately ``1e-16``, and a transformation
% back to ``0`` is done on retrieval.
%
% A square in the integer space is defined using three integers: the lowest x and y coordinates of the square, and
% the side length. This is the format of the squares to explore stored in the array `.todo`.
%
% In the sparse matrix `.checked`, we stored whether the square [x..x+s] x [y..y+s] has been checked by storing the
% side length ``s`` at the coordinate (x+1,y+1) (conversion from 0-based to 1-based).
%
% The class is constructed by asking
%

    properties (SetAccess = protected) % Immutable properties
        f % (function_handle): Function of two arguments to plot
        div % (integer(1,2)): Number of divisions on the x and y axis, must be a power of two
        axis % (double(1,4)): [xmin, xmax, ymin, ymax]
        positiveStyle % cell: If non-empty, plot positive points using the given plot arguments
        negativeStyle % cell: If non-empty, plot negative points using the given plot arguments
        boundaryStyle % cell: If non-empty, plot the boundary using the given plot arguments
    end

    properties (SetAccess = protected) % Mutable properties
        values % (double sparse matrix): Matrix of known function values, with evaluated zeros substituted by ``eps``
        checked % (integer sparse matrix)
        todo % (integer(3,\*)): Squares to explore, triplet (ix, iy, is) describing the square [ix,ix+is]*[iy,iy+is]
    end

    methods

        function self = ContourPlot(f, axis, pt, positiveStyle, negativeStyle, boundaryStyle)
        % Constructs an empty contour plot
        %
        % The contour plot is constructed using default div/step values. It adds in the ``todo`` squares list
        % the squares around the given point.
        %
        % Args:
        %   f (function_handle): Function taking two arguments (x,y) and returning a double value
        %   axis (double(1,4)): Values of [xmin, xmax, ymin, ymax], as with the ``axis`` Matlab command
        %   pt (double(1,2)): Point on the boundary, given using the double coordinates
            if nargin < 6
                boundaryStyle = {'k-'};
            end
            if nargin < 5
                negativeStyle = [];
            end
            if nargin < 4
                positiveStyle = [];
            end
            div = [2^18 2^18];
            s = 2^8;
            xmin = axis(1);
            xmax = axis(2);
            ymin = axis(3);
            ymax = axis(4);
            ix = (pt(1)-xmin)/(xmax-xmin)*div(1);
            iy = (pt(2)-ymin)/(ymax-ymin)*div(2);
            ix = round(ix/s)*s;
            iy = round(iy/s)*s;
            todo = [ix ix+s ix-s ix   ix
                    iy iy   iy   iy+s iy-s
                    s  s    s    s    s];
            self.f = f;
            self.div = div;
            self.axis = axis;
            self.todo = todo;
            self.values = sparse(div(1)+1, div(2)+1);
            self.checked = sparse(div(1)+1, div(2)+1);
            self.positiveStyle = positiveStyle;
            self.negativeStyle = negativeStyle;
            self.boundaryStyle = boundaryStyle;
        end

        function n = nx(self)
            n = self.div(1);
        end

        function n = ny(self)
            n = self.div(2);
        end

        function x = xdelta(self)
            x = self.axis(2) - self.axis(1);
        end

        function y = ydelta(self)
            y = self.axis(4) - self.axis(3);
        end

        function x = xstep(self)
            x = (self.axis(2) - self.axis(1))/self.div(1);
        end

        function y = ystep(self)
            y = (self.axis(4) - self.axis(3))/self.div(2);
        end

        function x = xmin(self)
            x = self.axis(1);
        end

        function x = xmax(self)
            x = self.axis(2);
        end

        function y = ymin(self)
            y = self.axis(3);
        end

        function y = ymax(self)
            y = self.axis(4);
        end

        function [z, flag] = evaluate(self, ix, iy)
        % Evaluates the function at the given integer index
        %
        % Returns an empty array if outside bounds
        %
        % Args:
        %   ix (integer): Integer zero-based x index on the grid
        %   iy (integer): Integer zero-based y index on the grid
        %
        % Returns
        % -------
        %   z: double
        %     Function value
        %   flag: integer
        %     Already known (0), just computed (1), outside bounds (-1)
            z = inf;
            if ix < 0 || ix > self.nx
                flag = -1;
                return
            elseif ix == 0
                x = self.xmin;
            elseif ix == self.nx
                x = self.xmax;
            else
                x = self.xmin + ix*self.xstep; % formula for the computation of x given ix
            end
            if iy < 0 || iy > self.ny
                flag = -1;
                return
            elseif iy == 0
                y = self.ymin;
            elseif iy == self.ny
                y = self.ymax;
            else
                y = self.ymin + iy*self.ystep; % formula for the computation of y given iy
            end
            z = full(self.values(ix+1, iy+1));
            if z ~= 0
                if z == eps
                    z == 0;
                end
                flag = 0;
            else
                flag = 1;
                f = self.f;
                z = f(x, y);
                %fprintf('Evaluting (%d,%d)=(%f,%f) => %f\n', ix, iy, x, y, z);
                if z == 0
                    self.values(ix+1, iy+1) = eps;
                else
                    self.values(ix+1, iy+1) = z;
                end
                if z >= 0
                    if ~isempty(self.positiveStyle)
                        plot(x, y, self.positiveStyle{:});
                    end
                else
                    if ~isempty(self.negativeStyle)
                        plot(x, y, 'x', self.negativeStyle{:});
                    end
                end
            end
        end


        function step(self)
        % Checks one squares of the todo list
            if isempty(self.todo)
                return
            end
            t = self.todo(:,end);
            self.todo = self.todo(:,1:end-1);
            s = t(3);
            ix1 = t(1);
            iy1 = t(2);
            if self.checked(ix1+1, iy1+1) == s
                % already checked
                return
            else
                self.checked(ix1+1, iy1+1) = s;
            end
            ix2 = t(1)+s;
            iy2 = t(2)+s;
            [z11, f11] = self.evaluate(ix1, iy1);
            [z12, f12] = self.evaluate(ix1, iy2);
            [z21, f21] = self.evaluate(ix2, iy1);
            [z22, f22] = self.evaluate(ix2, iy2);
            z = [z11 z12 z21 z22];
            f = [f11 f12 f21 f22];
            if any(f == -1)
                % outside bounds
                return
            end
            if all(z > 0) || all(z < 0)
                % all inside or all outside
                return
            end
            pts = [];
            if sign(z11) ~= sign(z12)
                pts = [pts; 0 -z11/(z12-z11)];
                self.todo = [self.todo [ix1-s; iy1; s]];
            end
            if sign(z12) ~= sign(z22)
                pts = [pts; -z12/(z22-z12) 1];
                self.todo = [self.todo [ix1; iy1+s; s]];
            end
            if sign(z21) ~= sign(z22)
                pts = [pts; 1 -z21/(z22-z21)];
                self.todo = [self.todo [ix1+s; iy1; s]];
            end
            if sign(z11) ~= sign(z21)
                pts = [pts; -z11/(z21-z11) 0];
                self.todo = [self.todo [ix1; iy1-s; s]];
            end
            if size(pts, 2) == 2 && ~isempty(self.boundaryStyle)
                X = self.xmin + (ix1 + pts(:,1) * s) * self.xstep;
                Y = self.ymin + (iy1 + pts(:,2) * s) * self.ystep;
                plot(X,Y,self.boundaryStyle{:});
            end
        end

    end

end
