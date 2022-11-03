clc,clear
syms sintheta3 costheta3 u sintheta2 costheta2 u2 sintheta1 costheta1 u1
xyzmatrix =  [0 0 1 -206;0 -1 0 0;1 0 0 0;0 0 0 1]*[0.5 0 -0.866 330;0 1 0 472;0.866 0 0.5 740-373;0 0 0 1];

that = [0 0 1 -206;0 -1 0 0;1 0 0 0;0 0 0 1];
trans = [0.5 0 -0.866 330;0 1 0 472;0.866 0 0.5 740-373;0 0 0 1];
thing = trans*that;
xyz = thing(1:3,4);
x = xyz(1);
y = xyz(2);
z = xyz(3);

%dhmatrix stdDH
[alpha0,a0,d1] = deal(0,0,0);
[alpha1,a1,d2] = deal(-90,-30,0);
[alpha2,a2,d3] = deal(0,340,0);
[alpha3,a3,d4] = deal(-90,-40,338);
[alpha4,a4,d5] = deal(90,0,0);
[alpha5,a5,d6] = deal(-90,0,0);

f1 = a3*costheta3 + d4*sind(alpha3)*sintheta3 +a2;
f2 = a3*cosd(alpha2)*sintheta3- d4*sind(alpha3)*cosd(alpha2)*costheta3-d4*sind(alpha2)*cosd(alpha3) - d3*sind(alpha2);
f3 = a3*sind(alpha2)*sintheta3 -d4*sind(alpha3)*sind(alpha2)*costheta3 +d4*cosd(alpha2) * cosd(alpha3) + d3 *cosd(alpha2);
    
k1 = f1;
k2 = -f2;
k3 = f1.^2 + f2.^2 + f3.^2+ a1.^2 + d2.^2 + 2*d2*f3;
k4 = f3*cosd(alpha1) + d2*cosd(alpha1);
r = x.^2 + y.^2 + z.^2;


%solve theta3 use r generate theta3
%enq1 = ((r-k3)^2)/(4*(-a1)^2) + ((z-k4)^2)/(sind(alpha1))^2 -k1^2 -k2^2;
enq = ((r-k3).^2)/(4*(a1).^2) + ((z-k4).^2)/(sind(alpha1)).^2 -k1.^2 -k2.^2 ==0;

replace_sin =  2*u/(1+u.^2);
replace_cos = (1-u.^2)/(1+u.^2);
enq_new1 = subs(enq,costheta3,replace_cos);
enq_new = subs(enq_new1,sintheta3,replace_sin);


result = solve(enq_new,u);

root = vpa(result,20);

angle = 2*atand(root);

angle_theta3 = vpa(angle,6);


%theta2
t3 = angle_theta3(3);
known_cost3 = cosd(t3);
known_sint3 = sind(t3);

f1_t2 = subs(f1, costheta3, known_cost3);
f1_t2 = subs(f1_t2, sintheta3, known_sint3);

f2_t2 = subs(f2, costheta3, known_cost3);
f2_t2 = subs(f2_t2, sintheta3, known_sint3);

f3_t2 = subs(f3, costheta3, known_cost3);
f3_t2 = subs(f3_t2, sintheta3, known_sint3);
    
k1_t2 = subs(k1, f1, f1_t2);
k2_t2 = subs(k2, f2, f2_t2);

k3_t2 = subs(k3, f1, f1_t2);
k3_t2 = subs(k3_t2, f2, f2_t2);
k3_t2 = subs(k3_t2, f3, f3_t2);

k4_t2 = subs(k4, f3, f3_t2);

replace2_sin =  2*u2/(1+(u2).^2);
replace2_cos = (1-(u2).^2)/(1+(u2).^2);

eqn_t2 = (k1_t2*costheta2 + k2_t2*sintheta2)*2*a1 + k3_t2 - r == 0;
eqn_t2 = subs(eqn_t2, costheta2, replace2_cos);
eqn_t2 = subs(eqn_t2, sintheta2, replace2_sin);

result_t2 = solve(eqn_t2,u2);
root_t2 = vpa(result_t2,20);
angle_t2 = 2*atand(root_t2);
angle_theta2 = vpa(angle_t2,6);

%theta1
t2 = angle_theta2(1);

g1 = cosd(t2)*f1_t2 - sind(t2)*f2_t2 + a1;
g2 = sind(t2)*cosd(alpha1)*f1_t2 + cosd(t2)*cosd(alpha1)*f2_t2 - sind(alpha1)*f3_t2 - d2*sind(alpha1);
g3 = sind(t2)*sind(alpha1)*f1_t2 + cosd(t2)*sind(alpha1)*f2_t2 + cosd(alpha1)*f3_t2 + d2*cosd(alpha1);

eqn_t1 = costheta1*g1 - sintheta1*g2 - x == 0;

replace1_sin =  2*u1/(1+(u1).^2);
replace1_cos = (1-(u1).^2)/(1+(u1).^2);

eqn_t1 = subs(eqn_t1, costheta1, replace1_cos);
eqn_t1 = subs(eqn_t1, sintheta1, replace1_sin);

result_t1 = solve(eqn_t1,u1);
root_t1 = vpa(result_t1,20);
angle_t1 = 2*atand(root_t1);
angle_theta1 = vpa(angle_t1,6);

t1 = angle_theta1(1);

%theta4 theta5 theta6
R01 = getR(0,0,0,double(t1));
R12 = getR(-90, -30, 0, double(t2));
R23 = getR(0, 340, 0, double(t3));

R03 = R01*R12*R23;

R06 = thing(1:3,1:3);
R03_p = R03*rotate_along_X(alpha3);
R36_p = R03_p\R06;



beta = atan2d( sqrt((R36_p(3,1))^2+(R36_p(3,2))^2) , R36_p(3,3));
alpha = atan2d( R36_p(2,3)/sind(double(beta)), R36_p(1,3)/sind(double(beta)));
gamma = atan2d( R36_p(3,2)/sind(double(beta)), -R36_p(3,1)/sind(double(beta)));

t4 = alpha + 180;
t5 = beta;
t6 = gamma + 180;

%t4 = alpha;
%t5 = beta;
%t6 = gamma;

R34 = getR(-90,-40,338,double(t4));
R45 = getR(90, 0, 0, double(t5));
R56 = getR(-90, 0, 0, double(t6));

R06_check = R01*R12*R23*R34*R45*R56;

function y = getR(alpha, a, d, theta)
    y = [cosd(theta) -sind(theta) 0; sind(theta)*cosd(alpha) cosd(theta)*cosd(alpha) -sind(alpha); sind(theta)*sind(alpha) cosd(theta)*sind(alpha) cosd(alpha) ];
end

function y = rotate_along_X(theta)
    y = [1 0 0; 0 cosd(theta) -sind(theta); 0 sind(theta) cosd(theta)];
end


