function img_fr = guide_2D_fr(img, points1, points2, faces)
points1 = points1';
points2 = points2';
img_fr = warping_image(img, points1, points2, faces');
end