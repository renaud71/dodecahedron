%maps spherical images on dodecahedron.
%prog., R. Toussaint, IPGS-EOST, CNRS Unistra, 2018.
clear all;
close all;

%colorcarte=imread('imagesourcered2.jpg'); % opening the image where color is stored as function ot longitude, latitude (image with horizontal pixels twice larger than vertical pixel size, 360 degrees versus 180)
colorcarte2=imread('export_fig_out_l5.png'); %data, C. Zaroli, Geoph. Journal International, 2016, DOI: 10.1093/gji/ggw315
colorcarte=colorcarte2(:,1:end-1,:); % particularity of tomo map with a 0 last column to remove
figure;
imshow(colorcarte);
pause(1);
gr=(1+sqrt(5))/2;  % golden ratio - nombre d or
vertex(1,:)=[1 1 1]; % coordinates of dodecahedron vertices - coordonnees de sommets du dodecahedre % see e.g. https://en.wikipedia.org/wiki/Regular_dodecahedron
vertex(2,:)=[1 1 -1];
vertex(3,:)=[1 -1 1];
vertex(4,:)=[1 -1 -1];
vertex(5,:)=[-1 1 1];
vertex(6,:)=[-1 1 -1];
vertex(7,:)=[-1 -1 1];
vertex(8,:)=[-1 -1 -1];
vertex(9,:)=[0 gr 1/gr];
vertex(10,:)=[0 gr -1/gr];
vertex(11,:)=[0 -gr 1/gr];
vertex(12,:)=[0 -gr -1/gr];
vertex(13,:)=[1/gr 0 gr];
vertex(14,:)=[1/gr 0 -gr];
vertex(15,:)=[-1/gr 0 gr];
vertex(16,:)=[-1/gr 0 -gr];
vertex(17,:)=[gr 1/gr 0];
vertex(18,:)=[gr -1/gr 0];
vertex(19,:)=[-gr 1/gr 0];
vertex(20,:)=[-gr -1/gr 0];
vertex=vertex./sqrt(dot(vertex(1,:),vertex(1,:))); % normalisation des vecteurs du centre aux sommets

%rotation pour mettre face sup à plat
vertexbis(:,2)=vertex(:,2).*0.85065080835204 + vertex(:,3).*0.5257311121191335527;
vertexbis(:,3)=-vertex(:,2).*0.5257311121191335527+vertex(:,3).*0.85065080835204;
vertex(:,2)=vertexbis(:,2);
vertex(:,3)=vertexbis(:,3);

X=vertex(:,1); % vecteurs avec les coords X, y, z des douze sommets
Y=vertex(:,2);
Z=vertex(:,3);
% scatter3(X,Y,Z); % representation 3D des vecteurs
% hold on
F(1,:)=[15 13 3 11 7];  % definition des faces comme suite de 5 sommets par leurs index, dans le sens des aiguilles d'une montre vu de l'exterieur - definition of faces, index list of vertices, counted clockwise seen from outside
F(2,:)=[13 1 17 18 3];
F(3,:)=[13 15 5 9 1];
F(4,:)=[5 15 7 20 19];
F(5,:)=[20 7 11 12 8];
F(6,:)=[12 11 3 18 4];
F(7,:)=[4 18 17 2 14];
F(8,:)=[2 17 1 9 10];
F(9,:)=[10 9 5 19 6];
F(10,:)=[6 19 20 8 16];
F(11,:)=[16 8 12 4 14];
F(12,:)=[16 14 2 10 6];

imax=12; % nb de faces - number of faces

for i=1:imax;
    % defining for each segment a vector perpendicular to the face from the
    % center to this segment, pointing outside the face - as a cross
    % product between two subsequent vertex vectors along the face edges
    for j=1:5
        if j<5
            n(i,j,:)=cross(vertex(F(i,j),:),vertex(F(i,j+1),:));
        else
            n(i,j,:)=cross(vertex(F(i,j),:),vertex(F(i,1),:));
        end
    end
    a(i,:)=vertex(F(i,5),:)-vertex(F(i,1),:); % a, vector from first vertex to last one in each face
    b(i,:)=vertex(F(i,2),:)-vertex(F(i,1),:); % b, vector from first vertex to second one in each face
    ua(i,:)=a(i,:)/sqrt(dot(a(i,:),a(i,:))); % ua, unit vector
    ub(i,:)=b(i,:)/sqrt(dot(b(i,:),b(i,:))); % unit vector
    c(i,:)=b(i,:)-(dot(ua(i,:),b(i,:))*ua(i,:)); % c, vector perp to a, inside the face
    uc(i,:)=c(i,:)/sqrt(dot(c(i,:),c(i,:))); % uc unit vector perp to ua
    uperp(i,:)=cross(ua(i,:),ub(i,:))/sqrt(dot(cross(ua(i,:),ub(i,:)),cross(ua(i,:),ub(i,:)))); % uperp, perp to ua,uc, direct orthonormal basis
end

lengthedge=sqrt(dot(a(1,:),a(1,:)));

Sizeimfinalx=2100;
Sizeimfinaly=2970;
L=250;
ca0=Sizeimfinalx*11/16;
cb0=Sizeimfinaly*4/16;

 figure;
 hold on;
for iface=1:12;
    if iface==1
        ori(iface,:)=[ca0 cb0];
        va(iface,:)=[-1 0]*L;
    end
    if iface==2
        ori(iface,:)=ori(iface-1,:)+vb(iface-1,:);
        va(iface,:)=ve(iface-1,:);
    end
    if iface==3
        ori(iface,:)=ori(iface-1,:);
        va(iface,:)=vb(iface-1,:);
    end
    if iface==4
        ori(iface,:)=ori(iface-1,:)+vb(iface-1,:)+ve(iface-1,:);
        va(iface,:)=vb(iface-1,:);
    end
    if iface==5
        ori(iface,:)=ori(iface-1,:)+vb(iface-1,:)+ve(iface-1,:)+vf(iface-1,:);
        va(iface,:)=ve(iface-1,:);
    end
    if iface==6
        ori(iface,:)=ori(iface-1,:)+vb(iface-1,:)+ve(iface-1,:)+vf(iface-1,:);
        va(iface,:)=ve(iface-1,:);
    end
     if iface==7
        ori(iface,:)=ori(iface-1,:)+va(iface-1,:);
        va(iface,:)=vf(iface-1,:);
     end
     if iface==8
         ori(iface,:)=ori(iface-1,:)+vb(iface-1,:)+ve(iface-1,:)+vf(iface-1,:);
         va(iface,:)=ve(iface-1,:);
     end
     if iface==9
         ori(iface,:)=ori(iface-1,:)+va(iface-1,:);
         va(iface,:)=vf(iface-1,:);
     end
     if iface==10
         ori(iface,:)=ori(iface-1,:)+va(iface-1,:);
         va(iface,:)=vf(iface-1,:);
     end
     if iface==11
         ori(iface,:)=ori(iface-1,:)+va(iface-1,:);
         va(iface,:)=vf(iface-1,:);
     end
     if iface==12
         ori(iface,:)=ori(iface-1,:);
         va(iface,:)=vg(iface-1,:);
     end
     
    vc(iface,:)=[-va(iface,2) va(iface,1)];
    vb(iface,:)=[va(iface,1)*cos(3*pi/5)-va(iface,2)*sin(3*pi/5) va(iface,2)*cos(3*pi/5)+va(iface,1)*sin(3*pi/5)];
    ve(iface,:)=-[vb(iface,1)*cos(3*pi/5)-vb(iface,2)*sin(3*pi/5) vb(iface,2)*cos(3*pi/5)+vb(iface,1)*sin(3*pi/5)];
    vf(iface,:)=-[ve(iface,1)*cos(3*pi/5)-ve(iface,2)*sin(3*pi/5) ve(iface,2)*cos(3*pi/5)+ve(iface,1)*sin(3*pi/5)];
    vg(iface,:)=-[vf(iface,1)*cos(3*pi/5)-vf(iface,2)*sin(3*pi/5) vf(iface,2)*cos(3*pi/5)+vf(iface,1)*sin(3*pi/5)];
    
    xa(iface,:)=[ori(iface,1) ori(iface,1)+va(iface,1)];
    ya(iface,:)=[ori(iface,2) ori(iface,2)+va(iface,2)];
    xb(iface,:)=[ori(iface,1) ori(iface,1)+vb(iface,1)];
    yb(iface,:)=[ori(iface,2) ori(iface,2)+vb(iface,2)];
    xe(iface,:)=[xb(iface,2) xb(iface,2)+ve(iface,1)];
    ye(iface,:)=[yb(iface,2) yb(iface,2)+ve(iface,2)];
    xf(iface,:)=[xe(iface,2) xe(iface,2)+vf(iface,1)];
    yf(iface,:)=[ye(iface,2) ye(iface,2)+vf(iface,2)];
    xg(iface,:)=[xf(iface,2) xf(iface,2)+vg(iface,1)];
    yg(iface,:)=[yf(iface,2) yf(iface,2)+vg(iface,2)];
    if iface==10
        xalang(iface,:)=[ori(iface,1) ori(iface,1)-(vb(iface,1)/5) ori(iface,1)+va(iface,1) ori(iface,1)+va(iface,1)];
        yalang(iface,:)=[ori(iface,2) ori(iface,2)-(vb(iface,2)/5) ori(iface,2)+va(iface,2) ori(iface,2)+va(iface,2)];
    else
        xalang(iface,:)=[ori(iface,1) ori(iface,1)-(vb(iface,1)/5) ori(iface,1)+va(iface,1)+(vg(iface,1)/5) ori(iface,1)+va(iface,1)];
        yalang(iface,:)=[ori(iface,2) ori(iface,2)-(vb(iface,2)/5) ori(iface,2)+va(iface,2)+(vg(iface,2)/5) ori(iface,2)+va(iface,2)];
    end
    if iface==1
        xblang(iface,:)=[ori(iface,1) ori(iface,1)-(va(iface,1)/5) ori(iface,1)+vb(iface,1) ori(iface,1)+vb(iface,1)];
        yblang(iface,:)=[ori(iface,2) ori(iface,2)-(va(iface,2)/5) ori(iface,2)+vb(iface,2) ori(iface,2)+vb(iface,2)];
    else
        xblang(iface,:)=[ori(iface,1) ori(iface,1)-(va(iface,1)/5) ori(iface,1)+vb(iface,1)-(ve(iface,1)/5) ori(iface,1)+vb(iface,1)];
        yblang(iface,:)=[ori(iface,2) ori(iface,2)-(va(iface,2)/5) ori(iface,2)+vb(iface,2)-(ve(iface,2)/5) ori(iface,2)+vb(iface,2)];
    end
    xelang(iface,:)=[xb(iface,2) xb(iface,2)+(vb(iface,1)/5) xb(iface,2)+ve(iface,1)-(vf(iface,1)/5)  xb(iface,2)+ve(iface,1)];
    yelang(iface,:)=[yb(iface,2) yb(iface,2)+(vb(iface,2)/5) yb(iface,2)+ve(iface,2)-(vf(iface,2)/5)  yb(iface,2)+ve(iface,2)];
    
    xflang(iface,:)=[xe(iface,2) xe(iface,2) xe(iface,2)+vf(iface,1)-(vg(iface,1)/5) xe(iface,2)+vf(iface,1)];
    yflang(iface,:)=[ye(iface,2) ye(iface,2) ye(iface,2)+vf(iface,2)-(vg(iface,2)/5) ye(iface,2)+vf(iface,2)];
    xflang(iface,:)=[xe(iface,2) xe(iface,2)+(ve(iface,1)/5) xe(iface,2)+vf(iface,1)-(vg(iface,1)/5) xe(iface,2)+vf(iface,1)];
    yflang(iface,:)=[ye(iface,2) ye(iface,2)+(ve(iface,2)/5) ye(iface,2)+vf(iface,2)-(vg(iface,2)/5) ye(iface,2)+vf(iface,2)];
    xglang(iface,:)=[xf(iface,2) xf(iface,2)+(vf(iface,1)/5) xf(iface,2)+vg(iface,1)+(va(iface,1)/5) xf(iface,2)+vg(iface,1)];
    yglang(iface,:)=[yf(iface,2) yf(iface,2)+(vf(iface,2)/5) yf(iface,2)+vg(iface,2)+(va(iface,2)/5) yf(iface,2)+vg(iface,2)];
    
    languetteindex(1,:)=[1 2 2 4 5];
    languetteindex(2,:)=[3 4 5 5 5];
    languetteindex(3,:)=[4 4 5 5 5];
    languetteindex(4,:)=[1 1 5 5 5];
    languetteindex(5,:)=[1 1 1 5 5];
    languetteindex(6,:)=[1 1 1 1 1];
    languetteindex(7,:)=[1 1 5 5 5];
    languetteindex(8,:)=[1 1 1 1 1];
    languetteindex(9,:)=[1 1 1 1 1];
    languetteindex(10,:)=[1 1 1 1 1];
    languetteindex(11,:)=[0 0 0 0 0];
    languetteindex(12,:)=[0 0 0 0 0];
    
    plot(xa(iface,:),ya(iface,:));
    plot(xb(iface,:),yb(iface,:));
    plot(xe(iface,:),ye(iface,:));
    plot(xf(iface,:),yf(iface,:));
    plot(xg(iface,:),yg(iface,:));
    
    if ismember(1,languetteindex(iface,:))
        plot(xalang(iface,:),yalang(iface,:));
    end
    if ismember(2,languetteindex(iface,:))
        plot(xblang(iface,:),yblang(iface,:));
    end
    if ismember(3,languetteindex(iface,:))
        plot(xelang(iface,:),yelang(iface,:));
    end
    if ismember(4,languetteindex(iface,:))
        plot(xflang(iface,:),yflang(iface,:));
    end
    if ismember(5,languetteindex(iface,:))
        plot(xglang(iface,:),yglang(iface,:));
    end
    %     plot(xa(iface,:),ya(iface,:));
    %     plot(xb(iface,:),yb(iface,:));
    %     plot(xe(iface,:),ye(iface,:));
    %     plot(xf(iface,:),yf(iface,:));
    %     plot(xg(iface,:),yg(iface,:));
    
end




Mmax=size(colorcarte,1);
Nmax=size(colorcarte,2);

% initializing a counter for all points counting, and one for counting the
% points in each face
count=0;
for i=1:imax
    countbis(i)=0;
end


for itheta=1:Mmax
    for iphi=1:Nmax
        count=count+1;
        phi=2*pi*(iphi)/Nmax;
        theta=pi*(itheta-1)/(Mmax-1);
        r=[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
        color=reshape(colorcarte(itheta,iphi,:),1,3);
        for i=1:imax;
            cond=1;
            for j=1:5
                if j<5
                    scal(i,j)=dot(reshape(n(i,j,:),1,3),r);
                    cond=cond*scal(i,j)<0;
                else
                    scal(i,j)=dot(reshape(n(i,j,:),1,3),r);
                    cond=cond*scal(i,j)<0;
                end
            end
            %rprojortho=r-(dot(r,uperp(i,:))*uperp(i,:));
            %rprojorthodec=rproj+(dot(vertex(F(i,1),:),uperp(i,:))*uperp(i,:));
            rprojdec=dot(vertex(F(i,1),:),uperp(i,:))/dot(r,uperp(i,:)).*r; %projection centrale sur la face
            if cond
                countbis(i)=countbis(i)+1;
                Xr(i,countbis(i))=rprojdec(1);
                Yr(i,countbis(i))=rprojdec(2);
                Zr(i,countbis(i))=rprojdec(3);
                Cr(i,countbis(i),:)=color(:);
                
%                 coordar=dot(rproj,ua(i,:));
%                 coordcr=dot(rproj,uc(i,:));
                coordar=dot(rprojdec,ua(i,:))-dot(vertex(F(i,1),:),ua(i,:));
                coordcr=dot(rprojdec,uc(i,:))-dot(vertex(F(i,1),:),uc(i,:));
                Planr(i,countbis(i),:)=ori(i,:) + coordar/lengthedge*va(i,:) + coordcr/lengthedge*vc(i,:);
            end
            Xu(count)=r(1);
            Yu(count)=r(2);
            Zu(count)=r(3);
        end
    end
end
% scatter3(Xu,Yu,Zu);
% hold on;
% scatter3(X,Y,Z);
figure;
hold on;
for i=1:imax;
    Couleur=reshape(Cr(i,:,:),size(Xr(i,:),2),3);
    Couleur=double(Couleur);
    Couleur=Couleur/255.;
    Xri=reshape(Xr(i,:),1,size(Xr,2));
    Yri=reshape(Yr(i,:),1,size(Xr,2));
    Zri=reshape(Zr(i,:),1,size(Xr,2));
    scatter3(Xri,Yri,Zri,80,Couleur,'filled');
end

%ExtPlanr=Planr;
%ExtPlanr(:,:,3)=0;
Xplanr(:,:)=Planr(:,:,1);
Yplanr(:,:)=Planr(:,:,2);

%Zplanr=0.*Yplanr;
% 
% figure;
% hold on;
% for i=1:imax;
%     Couleur2=reshape(Cr(i,:,:),size(Xplanr(i,:),2),3);
%     Couleur2=double(Couleur2);
%     Couleur2=Couleur2/255.;
%     Xplanri=reshape(Xplanr(i,:),1,size(Xplanr,2));
%     Yplanri=reshape(Yplanr(i,:),1,size(Xplanr,2));
%     Zplanri=reshape(Zplanr(i,:),1,size(Xplanr,2));
%     scatter3(Xplanri,Yplanri,Zplanri,80,Couleur2,'filled');
%     
% end


figure;
hold on;
for i=1:imax;
    Couleur2=reshape(Cr(i,:,:),size(Xplanr(i,:),2),3);
    Couleur2=double(Couleur2);
    Couleur2=Couleur2/255.;
    Xplanri=reshape(Xplanr(i,:),1,size(Xplanr,2));
    Yplanri=reshape(Yplanr(i,:),1,size(Xplanr,2));
    scatter(Xplanri,Yplanri,10,Couleur2,'filled');    
end
for iface=1:imax;
    plot(xa(iface,:),ya(iface,:));
    plot(xb(iface,:),yb(iface,:));
    plot(xe(iface,:),ye(iface,:));
    plot(xf(iface,:),yf(iface,:));
    plot(xg(iface,:),yg(iface,:));
    
        
    if ismember(1,languetteindex(iface,:))
        plot(xalang(iface,:),yalang(iface,:));
    end
    if ismember(2,languetteindex(iface,:))
        plot(xblang(iface,:),yblang(iface,:));
    end
    if ismember(3,languetteindex(iface,:))
        plot(xelang(iface,:),yelang(iface,:));
    end
    if ismember(4,languetteindex(iface,:))
        plot(xflang(iface,:),yflang(iface,:));
    end
    if ismember(5,languetteindex(iface,:))
        plot(xglang(iface,:),yglang(iface,:));
    end

end
axis equal;
axis off;
width=1500;
height=1500;
x0=0;
y0=0;
set(gcf,'units','points','position',[x0,y0,width,height]);
fig=gcf;
addpath ./altmany-export_fig-5b3965b/  %export matlab image package available at http://fr.mathworks.com/matlabcentral/fileexchange/23629-export_fig
export_fig(fig);


%% depth layers: 5,6,7,8,9,10,11,12: (a) 400?530; (b) 530?660; (c) 660?810; (d) 810?960; (e) 960?1110; (f) 1110?1310; (g) 1310?1510; (h) 1510?1710
%% central depth: 465, 595, 735, 885, 1035, 1210, 1410, 1610
%% radii (from 6371 km surface radius): 5906 km, 5776, 5636, 5486, 5336, 5161, 4961, 4761
