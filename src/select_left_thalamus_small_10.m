% Author: Veronica Corona, vc324@cam.ac.uk

function [roi]=select_left_thalamus_small_10(img)

if size(img,3)== 80
roi=img(95:124,80:118,35:47); % 36:46
else 

roi=img(70:94,94:129,126:149);


end
