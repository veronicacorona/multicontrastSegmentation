% Author: Veronica Corona, vc324@cam.ac.uk

function [roi]=select_left_thalamus_small(img)

if size(img,3)== 80
roi=img(95:124,80:118,35:47); % 36:46
else 
%   roi=img(70:99,91:129,126:149); % 30,39
    roi=img(73:97,94:129,126:149); % 36,25

end
