clear all;
display 'creating data set'
base = 128;
len1 = 4*base;
len2 = 2*base;
len3 = 4*base;

mat = 105*ones(len1,len2,len3);
mat = uint8(mat);

radius = [0.5*base, 0.5*base];
center = [[0.25*len1, 0.5*len2]; [0.75*len1, 0.5*len2]];
frequency = [[1/50, 1/75]; [1/50, 1/75]];
amplitude = [[0.05*len1, 0*len2]; [0.05*len1, 0*len2]];
phase = [[1, 2]; [1, 2]];
numChannels = length(radius);

for i=1:len1
	for j=1:len2
		for k=1:len3
			for c=1:numChannels
				if (i - center(c,1) - amplitude(c,1)*(cos(k*frequency(c,1)+phase(c,1))))^2 + ...
					 (j - center(c,2) - amplitude(c,2)*(cos(k*frequency(c,2)+phase(c,2))))^2 <= radius(c)^2
					mat(i,j,k) = 220;
				end
			end
		end
	end
end

disp 'done, writing out'
figure;
imshow(squeeze(mat(:,100,:)))
figure;
imshow(squeeze(mat(:,:,250)))

dicomwrite(reshape(mat,len1,len2,1,len3),'channels_512x256.dcm')

%gipl_write_volume(mat)
