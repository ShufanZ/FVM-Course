%% =====================
load("data.mat")


%% =======================subsonic

figure
hold on
for i=1:5
    plot(nstep,subRes(:,i),"lineWidth",3)
end
ylim([1e-5,1e1])
legend("adi","lj","pgs","pj","gslr")
set(gca,'yscale','log')
set(gca,'Fontsize',15)
xlabel("N step")
ylabel("Res")

%% =======================transonic

figure
hold on 
plot(nstep2,tranRes(:,1),"lineWidth",3)
plot(nstep2,tranRes(:,2),"lineWidth",3)
set(gca,'yscale','log')
legend("M=0.735","M=0.908")
set(gca,'Fontsize',15)
xlabel("N step")
ylabel("Res")
%% =======cp M=0.735
figure

plot(x,-cptrans1,"-o","lineWidth",3)
set(gca,'Fontsize',15)
xlabel("x/c")
ylabel("-c_p")
title("M_{inf}=0.735")
%% ======cp M = 0.908
figure
plot(x,-cptrans2,"-o","lineWidth",3)
set(gca,'Fontsize',15)
xlabel("x/c")
ylabel("-c_p")
title("M_{inf}=0.908")
%% ========M contour M = 0.908
for i=1:51
    for j=1:51
        fu(i,j)=F((i-1)*51+j,3);
        fm(i,j)=F((i-1)*51+j,5);
        fv(i,j)=F((i-1)*51+j,4);
    end
end
for i=1:51
    x(i) = F(i,1);
end
for j=1:51
    y(j) = F(1+51*(j-1),2);
end

starty = 0.1:0.2:2;
startx = -0.5*ones(size(starty));

figure
L = [1,1.1,1.2,1.3];
contour(x,y,fm,L,'lineWidth',3);
xlim([-0.5,1.5])
ylim([0,2])
set(gca,'Fontsize',15)
hold on
streamline(x,y,fu,fv,startx,starty)
xlabel("x/c")
ylabel("y/c")
%% ==========M contour M =0.735
for i=1:51
    for j=1:51
        fu1(i,j)=F2((i-1)*51+j,3);
        fm1(i,j)=F2((i-1)*51+j,5);
        fv1(i,j)=F2((i-1)*51+j,4);
    end
end
for i=1:51
    x1(i) = F2(i,1);
end
for j=1:51
    y1(j) = F2(1+51*(j-1),2);
end

starty = 0.1:0.2:2;
startx = -0.5*ones(size(starty));

figure
L = [0.75,0.8];
contour(x1,y1,fm1,L,'lineWidth',3);
xlim([-0.5,1.5])
ylim([0,2])
set(gca,'Fontsize',15)
hold on
streamline(x1,y1,fu1,fv1,startx,starty)
xlabel("x/c")
ylabel("y/c")