%H3 COMP 558
SurfFeaturepoints
%--------------------------Q1--------------------------------------------
bit = 1; %BIT TO REPRESENT WHICH IMAGE SET WE ARE DEALING WITH (1 or 2)
I1 = imread("p11.jpg"); %left image
I2 = imread("p12.jpg"); %right image

imtool(I1);
imtool(I2);

%Q1: Fundamental matrix estimation for Image pair 1
%    left         right
% (120,400) . (22,380) yes
% (410,226) . (356,220) yes
% (326,228) . (274,224) yes
% (680,370) . (418,222) yes
% (722,448) . (466,344) yes
% (514,238) . (462,90) yes
% (414,198) . (360,56) yes
% (595,138) . (636,120) yes

if (bit == 1)

    %List of X points in IM1 and IM2
    X1 = {120 410 326 680 722 514 414 595};
    X2 = {22 356 274 418 466 462 360 636};
    xave1 = (120 + 410 + 326 + 680 +722 +514 +414 + 595)/8;
    xave2 = (22 + 356 + 274 + 418 + 466 + 462 + 360 + 636)/8;

    %List of Y points in IM1 and IM2
    Y1 = {400 226 228 370 448 238 198 138};
    Y2 = {380 220 224 222 344 90 56 120};
    yave1 = (400 + 226 + 228 + 370 + 448 +238 + 198 + 138)/8;
    yave2 = (380 + 220 + 224 + 222+ 344 + 90 + 56 + 120)/8;

    %Initializing standard deviation values
    sd1 = 0;
    sd2 = 0;

    %Computing sd1 and sd2
    for i = 1:8
        sd1 = sd1 + ((X1{i} - xave1)^2 + (Y1{i} - yave1)^2);
        sd2 = sd2 + ((X2{i} - xave2)^2 + (Y2{i} - yave2)^2);
    end
    
    sd1 = (1/16) * sd1;
    sd2 = (1/16) * sd2;
    sd1 =  sqrt(sd1);
    sd2 =  sqrt(sd2);

    %Storing values before normalization
    bX1 = X1;
    bX2 = X2;
    bY1 = Y1;
    bY2 = Y2;
   
    %S1) Data Normalization     
    for i = 1:8
        X1{i} = (X1{i} - (xave1))/sd1;
        X2{i} = (X2{i} - (xave2))/sd2;
        Y1{i} = (Y1{i} - (yave1))/sd1;
        Y2{i} = (Y2{i} - (yave2))/sd2;
    end

    %S2) Least Squares
    %Initializing matrix A
    A = zeros(9,9);
    for r = 1:8
        A(r,:)=[(X1{r}*X2{r}) (Y1{r}*X2{r}) X2{r} (X1{r}*Y2{r}) (Y1{r}*Y2{r}) Y2{r} X1{r} Y1{r} 1];
    end

    %Using SVD matrix decomposition 
    [M, N, B] = svd(A);
    vector = B(:,9);
    F1 = [vector(1) vector(2) vector(3); vector(4) vector(5) vector(6); vector(7) vector(8) vector(9)];
   
    %S3) Mapping these entries to the entries of F via a de-normalization step
    M1 = [1/sd1 0 -xave1/sd1; 0 1/sd1 -yave1/sd1; 0 0 1];
    M2 = [1/sd2 0 -xave2/sd2; 0 1/sd2 -yave2/sd2; 0 0 1];
    F1 = M2'* F1 * M1;
   
    %S4) Enforcing a rank = 2 constraint on F
    [U, S, V] = svd(F1);
    S(3,3)= 0;
    F1 = U*S*V';
    %Displaying rank of matrix F1
    disp("Q1. Rank of Generated Fundamental Matrix is:");
    disp(rank(F1));
    %%%%GENERATING ESTIMATE OF FUNDAMENTAL MATRIX TO COMPARE RESULTS%%%%%%%
    left = [bX1{1} bY1{1}; bX1{2} bY1{2}; bX1{3} bY1{3}; bX1{4} bY1{4}; bX1{5} bY1{5}; bX1{6} bY1{6}; bX1{7} bY1{7}; bX1{8} bY1{8}];
    right = [bX2{1} bY2{1}; bX2{2} bY2{2}; bX2{3} bY2{3}; bX2{4} bY2{4}; bX2{5} bY2{5}; bX2{6} bY2{6}; bX2{7} bY2{7}; bX2{8} bY2{8}];
    f = estimateFundamentalMatrix(left,right,'Method','Norm8Point');
    
    disp("Q1. Fundamental Matrix for Image Pair #1 is:");
    disp(F1);

 else
    %{
    %Q1: Fundamental matrix estimation for Image pair 2
   %Left         Right
    (298 234)    (286 236)
    (526 390)    (416 390)
    (568 330)    (520 330)
    (362 86)     (350 90)
    (340 480)    (186 474)
    (78 378)     (48 376)
    (132 456)    (22 448)
    (430 330)    (406 330)
    %}

    %List of X points in IM1 and IM2
    X1 = {298 526 568 362 340 78 132 430};
    X2 = {286 416 520 350 186 48 22 406};
    xave1 = (298 + 526 + 568 + 362 + 340 + 78 + 132 +430)/8;
    xave2 = (286 + 416 + 520 + 350 + 186 + 48 + 22 + 406)/8;

    %List of Y points in IM1 and IM2
    Y1 = {234 390 330 86 480 378 456 330};
    Y2 = {236 390 330 90 474 376 448 330};
    yave1 = (234 +390+ 330 +86 +480 + 378 + 456+ 330)/8;
    yave2 = (236 +390 +330 + 90 + 474 + 376 +448 +330)/8;

    %Initializing standard deviation values
    sd1 = 0;
    sd2 = 0;

    %Computing sd1 and sd2
    for i = 1:8
        sd1 = sd1 + ((X1{i} - xave1)^2 + (Y1{i} - yave1)^2);
        sd2 = sd2 + ((X2{i} - xave2)^2 + (Y2{i} - yave2)^2);
    end
    
    sd1 = (1/16) * sd1;
    sd2 = (1/16) * sd2;
    sd1 =  sqrt(sd1);
    sd2 =  sqrt(sd2);

    %Storing values before normalization
    bX1 = X1;
    bX2 = X2;
    bY1 = Y1;
    bY2 = Y2;
   
    %S1) Data Normalization     
   for i = 1:8
    X1{i} = (X1{i} - (xave1))/sd1;
    X2{i} = (X2{i} - (xave2))/sd2;
    Y1{i} = (Y1{i} - (yave1))/sd1;
    Y2{i} = (Y2{i} - (yave2))/sd2;
   end

    %S2) Least Squares
    %Initializing matrix A
    A = zeros(9,9);
    for r = 1:8
        A(r,:)=[(X1{r}*X2{r}) (Y1{r}*X2{r}) X2{r} (X1{r}*Y2{r}) (Y1{r}*Y2{r}) Y2{r} X1{r} Y1{r} 1];
    end

    %Using SVD matrix decomposition 
    [M, N, B] = svd(A);
    vector = B(:,9);
    F1 = [vector(1) vector(2) vector(3); vector(4) vector(5) vector(6); vector(7) vector(8) vector(9)];
   
    %S3) Mapping these entries to the entries of F via a de-normalization step
    M1 = [1/sd1 0 -xave1/sd1; 0 1/sd1 -yave1/sd1; 0 0 1];
    M2 = [1/sd2 0 -xave2/sd2; 0 1/sd2 -yave2/sd2; 0 0 1];
    F1 = M2'* F1 * M1;
   
    %S4) Enforcing a rank = 2 constraint on F
    [U, S, V] = svd(F1);
    S(3,3)= 0;
    F1 = U*S*V';
    %Displaying rank of matrix F1
    disp("Q1. Rank of Generated Fundamental Matrix #2 is:");
    disp(rank(F1));
    %%%%GENERATING ESTIMATE OF FUNDAMENTAL MATRIX TO COMPARE RESULTS%%%%%%%
    left = [bX1{1} bY1{1}; bX1{2} bY1{2}; bX1{3} bY1{3}; bX1{4} bY1{4}; bX1{5} bY1{5}; bX1{6} bY1{6}; bX1{7} bY1{7}; bX1{8} bY1{8}];
    right = [bX2{1} bY2{1}; bX2{2} bY2{2}; bX2{3} bY2{3}; bX2{4} bY2{4}; bX2{5} bY2{5}; bX2{6} bY2{6}; bX2{7} bY2{7}; bX2{8} bY2{8}];
    f = estimateFundamentalMatrix(left,right,'Method','Norm8Point');
    disp("Q1. Fundamental Matrix for Image Pair #2 is:");
    disp(F1);
end

%--------------------------Q2--------------------------------------------
    %Extracting location of Surf Points (left and right)
    candidate_left = matchedLeft.Location;
    candidate_right =  matchedRight.Location;

    %Extracting 120 random points to iterate with RANSAC
    red_left = zeros(40,2);
    red_right =  zeros(40,2);
    random_40 = randperm(length(candidate_right),40);
    random_40 = random_40';
   
    for i = 1:length(random_40)
        red_left(i,:)= candidate_left(random_40(i,1),:);
        red_right(i,:)= candidate_right(random_40(i,1),:);   
    end
   
    %Want to find matrix F such that x2^T*F*x1 ≈ 0,using RANSAC
    %Implementing RANSAC
    %Adjusting Parameters for RANSAC
    consensus_lim = 8;

    %STORING 40 POINTS IN ARRAY IN CASE WE NEED THEM LATER
    all_left = red_left;
    all_right = red_right;
    F_r = [];
    F_l = [];
    iter_out = 8000;
    iter_in = 0;
    F_final= zeros(3,3);
    inlier = 0;
    random_nums = {};

   while (iter_in < iter_out)    
        %S1) Picking 8 pixel pairs randomly from the list of candidate pair 
        left_8 = zeros(8,2);
        right_8 =  zeros(8,2);
        random_8 = randperm(length(random_40),8);
        random_8 = random_8';
        random_nums{end+1} = random_8;

        for i = 1:8
            left_8(i,:)= all_left(random_8(i,1),:); 
            right_8(i,:)= all_right(random_8(i,1),:);
        end
        g = random_8';
     
        %S2) Estimating F (iter)
        %List of X points in IM1 and IM2
        X1 = {};
        X2 = {};
        Y1 = {};
        Y2 = {};

        for i = 1:8
            X1 {end +1} = left_8(i,1);
            X2 {end +1} = right_8(i,1);
            Y1 {end +1} = left_8(i,2);
            Y2 {end +1} = right_8(i,2);
        end
        
        %Getting X1, X2, Y1 and Y2 averages
        xave1 = 0;
        xave2 = 0;
        yave1 = 0;
        yave2 = 0;

        for i = 1:8
            xave1 = xave1 + X1{i};
            xave2 = xave2 + X2{i};
            yave1 = yave1 + Y1{i};
            yave2 = yave2 + Y2{i};
        end

        xave1 = xave1/8;
        xave2 = xave2/8;
        yave1 = yave1/8;
        yave2 = yave2/8;

        %Initializing standard deviation values
        sd1 = 0;
        sd2 = 0;
        %Computing sd1 and sd2
        for i = 1:8
            sd1 = sd1 + ((X1{i} - xave1)^2 + (Y1{i} - yave1)^2);
            sd2 = sd2 + ((X2{i} - xave2)^2 + (Y2{i} - yave2)^2);
        end

        sd1 = (1/16) * sd1;
        sd2 = (1/16) * sd2;
        sd1 =  sqrt(sd1);
        sd2 =  sqrt(sd2);

        %Storing values before normalization
        bX1 = X1;
        bX2 = X2;
        bY1 = Y1;
        bY2 = Y2;
   
        %2.S2) Least Squares
        %Initializing matrix A
        A = zeros(9,9);
        for r = 1:8
            A(r,:)=[(X1{r}*X2{r}) (Y1{r}*X2{r}) X2{r} (X1{r}*Y2{r}) (Y1{r}*Y2{r}) Y2{r} X1{r} Y1{r} 1];
        end

        %Using SVD matrix decomposition 
        [M, N, B] = svd(A);
        vector = B(:,9);
        F = [vector(1) vector(2) vector(3); vector(4) vector(5) vector(6); vector(7) vector(8) vector(9)];    
   
        %2.S3) Mapping these entries to the entries of F via a de-normalization step
        M1 = [1/sd1 0 -xave1/sd1; 0 1/sd1 -yave1/sd1; 0 0 1];
        M2 = [1/sd2 0 -xave2/sd2; 0 1/sd2 -yave2/sd2; 0 0 1];
        F = M2'*F* M1;
   
        %2.S4) Enforcing a rank = 2 constraint on F
        [U, S, V] = svd(F);
        S(3,3)= 0;
        F = U*S*V';
       
    %S3) Determining which remaining points in set that belong to Consensus set
        cX1 = {};
        cX2 = {};
        cY1 = {};
        cY2 = {};


        for j = 1:length(red_left)
           P2 = [red_right(j,1);red_right(j,2);1];
           P1 = [red_left(j,1);red_left(j,2);1];
           T = P2'*F*P1;   
           T = abs(T);
           
           if T < 0.006
                cX1 {end+1} = red_left(j,1);
                cX2 {end+1} = red_right(j,1);
                cY1 {end+1} = red_left(j,2);
                cY2 {end+1} = red_right(j,2);
           end
            
        end

    %If your consensus set is above a certain threshold, we store model of
    %F
      if ( (length(cX1) >= consensus_lim ) && length(cX1)> inlier)        
            F_final = F;
            inlier = length(cX1);   
            F_l = left_8; 
            F_r = right_8; 
    %Storing consensus set
            i_X1 = cX1;
            i_X2 = cX2;
            i_Y1 = cY1;
            i_Y2 = cY2;
            
      end 
    iter_in = iter_in +1;
   end

%-------------------------------Q3-----------------------------------------
%Finding e1 (right null space) and e2 (left null space) of previously generated F2 with RANSAC.

if (bit == 1)

F_final = [-1.22385320147802e-10	-2.50616093221418e-08	0.000451599943821839;
7.81213585694305e-09	1.81365234376140e-08	0.00617765081138932;
-0.000203816506415516	-0.00776863296924383	0.291008688968557];

else

F_final = [-1.53462543900070e-09	-3.81367453594227e-08	0.000959839530888257;
4.32486782071492e-08	-7.73082718392819e-09	-0.00701186455415368;
-0.000779682823153024	0.00697016201977165	-0.0364890075064721];

end

e1 = null(F_final);
e2 = null(F_final');


%X values range from 1 to 752
%Y values range from 1 to 500

H1 =[1 0 0;-e1(2)/e1(1) 1 0;-1/e1(1) 0 1];
H2 =[1 0 0;-e2(2)/e2(1) 1 0;-1/e2(1) 0 1];

%Initializing array to store both rectified images
copy_final = zeros(1200,2000,3, 'uint8');

%Initializing array to store rectified image 1
I1_new = zeros(1200,1000,3, 'uint8');

%Initializing array to store rectified image 2
I2_new = zeros(1200,1000,3, 'uint8');

for i =1:500 %ITERATING THROUGH Y VALUES  
    for j=1:752 %ITERATING THROUGH X VALUES 

     pixel_value_I1_r = I1(i,j,1);
     pixel_value_I1_g = I1(i,j,2);
     pixel_value_I1_b = I1(i,j,3);


     pixel_value_I2_r = I2(i,j,1);  
     pixel_value_I2_g = I2(i,j,2); 
     pixel_value_I2_b = I2(i,j,3); 
    
     point_I1 = [j;i;1];  
     point_I2 = [j;i;1];  
     
     new_point_I1 = H1 * point_I1;
     new_point_I2 = H2 * point_I2;

     new_point_I1 = round(new_point_I1);
     new_point_I2 = round(new_point_I2);
    
      if (new_point_I1(1,1) >= -299 && new_point_I1 (1,1) <= 1700 && new_point_I1 (2,1) >= -299 && new_point_I1 (2,1) <= 1700)
        I1_new(new_point_I1(2,1)+300,new_point_I1(1,1)+300,:) =  [pixel_value_I1_r;pixel_value_I1_g;pixel_value_I1_b];
        copy_final(new_point_I1(2,1)+300,new_point_I1(1,1)+300,:) =  [pixel_value_I1_r;pixel_value_I1_g;pixel_value_I1_b];
     else

     end

     if (new_point_I2(1,1) >= -299 && new_point_I2 (1,1) <= 1700 && new_point_I2(2,1) >= -299 && new_point_I2(2,1) <= 1700)
        I2_new(new_point_I2(2,1)+300,new_point_I2(1,1)+300,:) =  [pixel_value_I2_r;pixel_value_I2_g;pixel_value_I2_b];
        copy_final(new_point_I2(2,1)+300,new_point_I2(1,1)+1100,:) =  [pixel_value_I2_r;pixel_value_I2_g;pixel_value_I2_b];
     else
       
     end     

    end
end
    

%Displaying results

figure;
imshow(copy_final);


%-------------------------------Q4-----------------------------------------

if (bit == 1)
%5 selected points to get reconstruction from Image pair 1
%{
%Left         Right
(332 216)    (280 212)
(514 348)    (466 344)
(306 108)    (462 98)
(194 206)    (92 204)
(478 470)    (202 464)
%}

lef = [332 216; 514 348; 306 108; 194 206; 478 470];
rig = [280 212; 466 344; 462 98; 92 204; 202 464];

lefx = {332,514,306,194,478};
lefy = {216, 348, 108, 206, 470};

rigx = {280, 466, 462, 92, 202};
rigy = {212, 344, 98, 204, 464};

F_final = [-1.22385320147802e-10	-2.50616093221418e-08	0.000451599943821839;
7.81213585694305e-09	1.81365234376140e-08	0.00617765081138932;
-0.000203816506415516	-0.00776863296924383	0.291008688968557];

f = 31;

else

 %5 selected points to get reconstruction from image pair 2
    %{
    %Left         Right
    (298 234)    (286 236)
    (526 390)    (416 390)
    (568 330)    (520 330)
    (362 86)     (350 90)
    (340 480)    (186 474)
    %}

    lef = [298 234; 526 390; 568 330; 362 86; 340 480];
    rig = [286 236; 416 390; 520 330; 350 90; 186 474];

    lefx = {298,526,568,362,340};
    lefy = {234,390,330,86,480};

    rigx = {280, 416, 520, 350, 186};
    rigy = {236, 390, 330, 90, 474};

    F_final = [-1.53462543900070e-09 -3.81367453594227e-08 0.000959839530888257;
    4.32486782071492e-08 -7.73082718392819e-09 -0.00701186455415368;
    -0.000779682823153024 0.00697016201977165 -0.0364890075064721];

    f = 32;

end


sx = 0.0312; % = 1/mx
sy = 0.0312; % = 1/my


%Constructing matrix K

K = [f/sx 0 376;0 f/sy 250;0 0 1];


E =K'*F_final*K;

%Correcting E

[U,S,V] = svd(E);
S = [1 0 0; 0 1 0; 0 0 0];
E = U*S*V;


%We define matrix W to denote a vector orthogonal to l and r
W=[0 -1 0;1 0 0;0 0 1];

%Extracting Translation vector as right null space of E
T1 = null(E);
T2 = -null(E);

%So now we have that R options are equal to:

R1 = U*W*V';
R2 = U*W'*V';

if(det(R1) == -1)
    R1 = -(R1);
end
if (det(R2) == -1)
    R2 = -(R2);
end

figure;
imshow(I1);

for i =1:5

    pl = [lefx{i}; lefy{i}; f];
    pr = [rigx{i}; rigy{i}; f];

    Rtpr1 = R1'*pr;
    Rtpr2 = R2'*pr;
    
    %generating w1 and w2
    w1 = cross(pl, Rtpr1);
    w2 = cross(pl, Rtpr2);
    %generating w1 and w2

    c1 = inv([pl -Rtpr1 w1]) * T1;
    c2 = inv([pl -Rtpr1 w2]) * T1;
    c3 = inv([pl -Rtpr2 w1]) * T1;
    c4 = inv([pl -Rtpr2 w2]) * T1;

    c1_2 = inv([pl -Rtpr1 w1]) * T2;
    c2_2 = inv([pl -Rtpr1 w2]) * T2;
    c3_2 = inv([pl -Rtpr2 w1]) * T2;
    c4_2 = inv([pl -Rtpr2 w2]) * T2;

    %Computing all possible mid points, and plotting best one to image 
    apl_1_1 = c1(1,1)*pl; 
    br_1_1 = T1 + c1(2,1)*Rtpr1;  
    midpoint1_1 = (br_1_1+ apl_1_1)/2;

    apl_2_1 = c2(1,1)*pl;
    br_2_1 = T1  + c2(2,1)*Rtpr1;
    midpoint2_1 = (br_2_1+ apl_2_1)/2;

    apl_3_1 = c3(1,1)*pl;
    br_3_1 = T1 + c3(2,1)*Rtpr2;
    midpoint3_1 = (br_3_1+ apl_3_1)/2;
    
    apl_4_1 = c4(1,1)*pl;
    br_4_1 = T1+ c4(2,1)*Rtpr2;
    midpoint4_1 = (br_4_1+ apl_4_1)/2;
   
    %Best results
    apl_1_2 = c1_2(1,1)*pl;
    br_1_2 = T2 + c1_2(2,1)*Rtpr1;
    midpoint1_2 = (br_1_2+ apl_1_2)/2;

    %Best results too
    apl_2_2 = c2_2(1,1)*pl;
    br_2_2 = T2+ c2_2(2,1)*Rtpr1;
    midpoint2_2 = (br_2_2+ apl_2_2)/2;
    
    apl_3_2 = c3_2(1,1)*pl;  
    br_3_2 = T2 + c3_2(2,1)*Rtpr2;
    midpoint3_2 = (br_3_2+ apl_3_2)/2;

    apl_4_2 = c4_2(1,1)*pl;
    br_4_2 = T2 + c4_2(2,1)*Rtpr2;
    midpoint4_2 = (br_4_2+ apl_4_2)/2;
    if (bit == 2)
    txt = sprintf('%f, %d', abs(midpoint3_1(3,1)), i);
    else
    txt = sprintf('%f, %d', abs(midpoint1_2(3,1)), i);
    end
    text(lefx{i}, lefy{i},txt);
end


%-------------------------------Q5-----------------------------------------
%If we are dealing with image pair 1 we set this to 1 or 2 if image pair 2
bit  = 1;
if (bit ==1)
    f = 31;
else
    f = 32;
end

%Establishing Lukas-Kanade parameters
Tx = 100; %Tx parameter
s = 40; %Shift parameter
window = 11;
window_h = ceil(window/2);
%Initializing empty disparity array
disparity = [];

%Initializing empty depth map
depth = [];

%Using rectified images
left = I1_new;
right = I2_new;

%Using a 1D Lucas-Kanade type method along row x.

for i = 1+window:1200-window %Iterating though y values                                             
    for j = 1+window:1052-window %Iterating though x values                                           
        
        % Initializing window
        b_y = i-window_h:i+window_h;
        b_x = j-window_h:j+window_h;
        
        % Left image intensity at current point
        intensity_l = left(b_y,b_x);
        
        if (j+s < 1052-window_h)
            %Initializing counter
            temp = 1;
            array = [];
            % Sliding window (Lukas-Kanade).
            for k= j:j+s
                % Updating window values
                new_window = k-window_h:k+window_h;
                % Right image intensity at current point
                intensity_r = right(b_y,new_window,1); 
                % Computing disparity for each slide and store it.         
                array(temp) = sum(abs(intensity_l-intensity_r),'all'); 
                temp = temp+1;
            end
        
        end
        
        %Now we look for the minimum element of the array and store the index
        [~,m] = min(array);  

        %Storing value in the disparity array
        disparity(i,j) = m/s;  
   
        % Storing value in the depth map
        depth(i,j) = round((f*Tx)/disparity(i,j));
        
    end
end

depth = uint16(depth);
                                    
disparity = imgaussfilt(disparity,2);
depth = imgaussfilt(depth,2);

figure;
imshow(disparity);  
figure;
imshow(depth);






