function rotation = ntcRotation(state)
%NTCROTATION  Build the radial-tangential-cross rotation matrix.
%
%   Rows are unit radial, tangential, and cross-track directions from r, v.
%
%   Author: Pasquale Marzaioli

normal = state(1:3) / norm(state(1:3));
crossTrack = cross(state(1:3), state(4:6));
crossTrack = crossTrack / norm(crossTrack);
tangential = cross(crossTrack, normal);
rotation = [normal, tangential, crossTrack].';
end
