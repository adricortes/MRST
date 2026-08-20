classdef delaunayTriangulation
% Minimal Octave replacement for MATLAB's delaunayTriangulation class.
% Supports only what upr/clippedPebi2D needs: construction from a 2D
% point set and the edges() method.
   properties
      Points
      ConnectivityList
   end
   methods
      function obj = delaunayTriangulation(varargin)
         if nargin == 1
            p = varargin{1};
         else
            p = [varargin{1}(:), varargin{2}(:)];
         end
         obj.Points = p;
         obj.ConnectivityList = delaunay(p(:,1), p(:,2));
      end
      function E = edges(obj)
         t = sort(obj.ConnectivityList, 2);
         E = unique([t(:,[1 2]); t(:,[2 3]); t(:,[1 3])], 'rows');
      end
   end
end
