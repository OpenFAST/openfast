function DirCos = DirCosMtrx0( x0, y0, z0, x1, y1, z1 )
%DIRCOSMTRX Summary of this function goes here
%   x is the local axial axis.

DirCos = zeros(3, 3);

xz = sqrt((x0-x1)^2+(z0-z1)^2);
xyz = sqrt((x0-x1)^2+(y0-y1)^2+(z0-z1)^2);
if( xz==0 )
    if (y1<y0)
        DirCos = [0 1 0;
                  1 0 0;
                  0 0 -1];
    else
        DirCos = [0 1 0;
                  -1 0 0;
                  0 0  1];
    end
else
    DirCos(1, 3) = (x0-x1)*(y0-y1)/(xz*xyz);
    DirCos(1, 2) = (z0-z1)/xz;
    DirCos(1, 1) = (x1-x0)/xyz;

    DirCos(2, 3) = -xz/xyz;
    DirCos(2, 1) = (y1-y0)/xyz;

    DirCos(3, 3) = (y0-y1)*(z0-z1)/(xz*xyz);
    DirCos(3, 2) = (x1-x0)/xz;
    DirCos(3, 1) = (z1-z0)/xyz;
    
    
    
end

end

