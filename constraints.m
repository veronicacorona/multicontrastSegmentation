%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

function result = constraints(projector, center, radius)
    result = [];
    result.project = projector;
    if (exist('center','var'))
        result.center = center;
    end
    if (exist('radius','var'))
        result.radius = radius;
    end
end
