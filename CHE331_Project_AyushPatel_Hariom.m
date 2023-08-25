Ca=[0.5 3 1 2 6 10 8 4];
Ca0=[2 5 6 6 11 14 16 24];
Tau=[30 1 50 8 4 20 20 4];
rA=[];
% CA_in=input('input concenteration : ');
% CA_out=input('output concenteration : ');
% flow_rate=input('flow rate : ');

for i=1:8
    rate=Tau(i)/(Ca0(i)-Ca(i));
    rA=[rA rate];
end
% Create spline interpolation object
spl = csape(Ca, rA, 'clamped');

% Define new x values to interpolate
x_interp = linspace(0.5, 10, 1000);

% Evaluate the spline at the new x values
y_interp = ppval(spl, x_interp);

% Plot the original data and the interpolated spline
plot(Ca, rA, 'o', x_interp, y_interp, '-');
grid on;
% rate is second order (-ra)=k*Ca*Ca


% Single PFR
F=@(Ca) ppval(spl, Ca);
% spl = csape(Ca, rA, 'clamped');
x=linspace(1,10,9001);
y = ppval(spl, x);
area = trapz(x, y);
% fprintf('Area under the curve under %.4f\n',area);
Volume_PFR = area*0.1;
fprintf('Volume_PFR %.4f\n',Volume_PFR);

% Single MFR
area_cstr = y(1)*(max(Ca)-1);
% fprintf('Area under the curve under %.4f\n',area_cstr);
Volume_CSTR = area_cstr*0.1;
fprintf('Volume_CSTR %.4f\n',Volume_CSTR);

% Two MFR
% F=@(Ca) ppval(spl, Ca);
optimized_area=1e5;
Ca_opt=min(Ca);
k=1;
for j=1:0.001:10
    mfr_area = (j-min(x))*F(1) + (max(x)-x(k))*F(j);  % current area under curve
    optimized_area=min(optimized_area, mfr_area);
    if optimized_area==mfr_area
        Ca_opt=j;
    end
    k=k+1;
end
fprintf('Volume_CSTR in series %.4f\n',optimized_area*0.1);

% PFR + MFR
Ca_min_val=fminbnd(F,Ca(1),Ca(end));
rA_min=F(Ca_min_val);
x_1=linspace(x(1),Ca_min_val,1000);
y_1=zeros(1,1000);
for i=1:1000
    y_1(i)=F(x_1(i));
end
x_2=linspace(Ca_min_val, max(Ca), 1000);
y_2=ones(1, 1000)*rA_min;
area_pfr_mfr=  integral(F,x(1),Ca_min_val) + ((max(Ca)-Ca_min_val)*rA_min); %trapz(x_1,y_1) 
fprintf('Volume_PFR+MFR %.4f\n',area_pfr_mfr*0.1);

% PFR + Recycle
MIN_Val=1e5;
Ca_val=1;
tolrance=0.1;
i=1;
t=1;
x_3=linspace(1,10,1000);
while(true)
    val_avg = integral(F,x_3(1),x_3(i))/((x_3(i)-x_3(1)));
    Ca_val=x_3(i);
    if abs(F(x_3(i))-val_avg)<tolrance
          break;
    end
    
    i=i+1;
end
recycle_ratio=(max(x)-Ca_val)/(Ca_val-min(x));
fprintf('Volume_PFR with recycle %.4f\n',F(x_3(i))*(max(x)-min(x))*0.1);
