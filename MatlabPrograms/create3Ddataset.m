len1 = 500;
len2 = 200;
len3 = 500;

mat = 105*ones(len1,len2,len3);
mat = uint8(mat);

%zf = @(x) 2 * ( cos(x/8) + 10);
zf1 = @(x) cos(x/50) + 2;
zf2 = @(x) cos(x/50) + 2;

R1 = 50;
R2 = 100;

xshift1 = 50;
yshift1 = 50;

xshift2 = 150;
yshift2 = 150;

for i=1:len1
	for j=1:len2
		for k = 1:len3
			if ( i - xshift1 * zf1(k) ) ^ 2 + (j - yshift1) ^ 2 <= R1 * R1
				mat(i,j,k) = 220;
			end
			if ( i - xshift2 * zf2(k) ) ^ 2 + (j - yshift2) ^ 2 <= R2 * R2
				mat(i,j,k) = 220;
			end
		end
	end
end

imshow(squeeze(mat(:,35,:)))
imshow(squeeze(mat(:,:,30)))

