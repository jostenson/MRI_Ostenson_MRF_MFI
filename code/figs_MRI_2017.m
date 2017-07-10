%% MRI journal figures

clear, close all

MFI_thresh = 3;

%% Figure 1 PIQT load


datadir = '../data_out/';
datadir_B0 = '../data_in/';
fig_out = '../figures/';

fn1 = 'MRF_proc_PIQT.mat';
fn1MFI = 'MRF_proc_PIQT_MFI.mat';
fn1B0 = 'MRF_B0_PIQT.mat';

% load
load([datadir fn1])

img_mag = abs(sum(input_MRF_dp.img_AR_comb,3));
effMtx = size(img_mag,1);

mythresh_v = multithresh(img_mag,MFI_thresh);
my_T1_map = output_MRF_match.T1_map;
my_T2_map = output_MRF_match.T2_map;
my_T1_map(img_mag < mythresh_v(1)) = 0;
my_T2_map(img_mag < mythresh_v(1)) = 0;
T1_map_masked = my_T1_map;
T2_map_masked = my_T2_map;

load([datadir fn1MFI])
load([datadir_B0 fn1B0])

img_mag_MFI = abs(sum(input_MRF_dp.img_AR_comb,3));
B0_map_masked = output_B0.B0_map_masked;
mythresh_v = multithresh(img_mag_MFI,MFI_thresh);
my_T1_map = output_MRF_match.T1_map;
my_T2_map = output_MRF_match.T2_map;
my_T1_map(img_mag < mythresh_v(1)) = 0;
my_T2_map(img_mag < mythresh_v(1)) = 0;
T1_map_masked_MFI = my_T1_map;
T2_map_masked_MFI = my_T2_map;


%% Figure 1 PIQT create figure

% desired image resolution in dpi
resImg = 500;

% font sizes
fontSzFigLet = 40; % figure letter size
fontSzCBlab = 28; % colorbar label
fontSzCBTicks = 28; % colorbar tick labels

figure(1); clf;
imagesc(B0_map_masked);
set(gcf,'Color','w');
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
c = colorbar('FontSize',fontSzCBTicks,'LineWidth',2.0);
c.Label.String = 'B0 off-resonance (Hz)';
c.Label.FontSize = fontSzCBlab;
c.Position(1) = 0.83;
c.Position(2) = 0.15;
c.Position(3) = 0.03;
c.Position(4) = 0.75;
text(25,25,'c','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgB0 = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc(img_mag);
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(25,25,'a','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgNoMFI = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(3); clf;
imagesc(img_mag_MFI);
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(25,25,'b','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgMFI = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

[rowB0,colB0]= size(imgB0);
[rowImg,colImg] = size(imgNoMFI);

imgFig1 = [imgNoMFI  imgMFI  imgB0];
figure(4); clf;
imshow(imgFig1)
txtexec = sprintf('export_fig([fig_out ''MRI_Figure1_MRFMFI.png''],''-r%d'')',resImg);
eval(txtexec);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure1_MRFMFI.tif''],''-r%d'')',resImg);
eval(txtexec);

%% Figure 2 50 mL phantom load data
close all, clc

% files to be loaded
fn1 = 'MRF_proc_shim_exp_ref.mat';
fn1MFI = 'MRF_proc_shim_exp_ref_MFI.mat';
fn1B0 = 'MRF_B0_shim_exp_ref.mat';

fn2 = 'MRF_proc_shim_exp_y.mat';
fn2MFI = 'MRF_proc_shim_exp_y_MFI.mat';
fn2B0 = 'MRF_B0_shim_exp_y.mat';

fn3 = 'MRF_proc_shim_exp_yz1.mat';
fn3MFI = 'MRF_proc_shim_exp_yz1_MFI.mat';
fn3B0 = 'MRF_B0_shim_exp_yz1.mat';

fn4 = 'MRF_proc_shim_exp_z.mat';
fn4MFI = 'MRF_proc_shim_exp_z_MFI.mat';
fn4B0 = 'MRF_B0_shim_exp_z.mat';

img_mag = zeros(effMtx,effMtx,4);
T1_map_masked = img_mag;
T2_map_masked = img_mag;
img_mag_MFI = img_mag;
T1_map_masked_MFI = img_mag;
T2_map_masked_MFI = img_mag;
B0_map_masked = img_mag;

for ii = 1:4
    figure(100+ii); clf;

    txt = sprintf('load([datadir fn%d])',ii);
    eval(txt);
    
    my_img = abs(sum(input_MRF_dp.img_AR_comb,3));
    img_mag(:,:,ii) = my_img;
    mythresh_v = multithresh(my_img,MFI_thresh);
    my_T1_map = output_MRF_match.T1_map;
    my_T2_map = output_MRF_match.T2_map;
    my_T1_map(my_img < mythresh_v(1)) = 0;
    my_T2_map(my_img < mythresh_v(1)) = 0;
    T1_map_masked(:,:,ii) = my_T1_map;
    T2_map_masked(:,:,ii) = my_T2_map;
    
    txt = sprintf('load([datadir fn%dMFI])',ii);
    eval(txt);
    txt = sprintf('load([datadir_B0 fn%dB0])',ii);
    eval(txt);
    
    my_img = abs(sum(input_MRF_dp.img_AR_comb,3));
    img_mag_MFI(:,:,ii) = my_img;
    mythresh_v = multithresh(my_img,MFI_thresh);
    my_T1_map = output_MRF_match.T1_map;
    my_T2_map = output_MRF_match.T2_map;
    my_T1_map(my_img < mythresh_v(1)) = 0;
    my_T2_map(my_img < mythresh_v(1)) = 0;
    T1_map_masked_MFI(:,:,ii) = my_T1_map;
    T2_map_masked_MFI(:,:,ii) = my_T2_map;
    B0_map_masked(:,:,ii) = output_B0.B0_map_masked;
      
end

%% Figure 2 50 mL phantom create figure

% desired image resolution in dpi
resImg = 500;

% font sizes
fontSzFigLet = 40; % figure letter size
fontSzCBlab = 28; % colorbar label
fontSzCBTicks = 28; % colorbar tick labels

% T1 and T2 display ranges
minT1 = 0;
maxT1 = 2500;
minT2 = 0;
maxT2 = 500;

% clip sides of images
c1 = 21;
c2 = 220;
B0_map_masked = B0_map_masked(:,c1:c2,:);
T1_map_masked = T1_map_masked(:,c1:c2,:);
T2_map_masked = T2_map_masked(:,c1:c2,:);
T1_map_masked_MFI = T1_map_masked_MFI(:,c1:c2,:);
T2_map_masked_MFI = T2_map_masked_MFI(:,c1:c2,:);

% row 1
r = 1;

figure(1); clf;
imagesc(B0_map_masked(:,:,r));
caxis([-450 450]);
set(gcf,'Color','w');
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
c = colorbar('north','FontSize',fontSzCBTicks,'LineWidth',2.0);
c.Label.String = 'B0 off-resonance (Hz)';
c.Label.FontSize = fontSzCBlab;
c.Color = 'w';
c.Position(1) = 0.30;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
text(25,190,'a','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('imgB0_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc(T1_map_masked(:,:,r)); caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
text(25,25,'w/o MFI','FontSize',2*fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T1NoMFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(3); clf;
imagesc(T1_map_masked_MFI(:,:,r)); caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
text(65,25,'MFI','FontSize',2*fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T1MFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(4); clf;
imagesc(T2_map_masked(:,:,r)); caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
text(25,25,'w/o MFI','FontSize',2*fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T2NoMFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(5); clf;
imagesc(T2_map_masked_MFI(:,:,r)); caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
text(65,25,'MFI','FontSize',2*fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T2MFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

% row 2
r = 2;

figure(1); clf;
imagesc(B0_map_masked(:,:,r));
caxis([-450 450])
set(gcf,'Color','w');
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(25,190,'b','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('imgB0_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc(T1_map_masked(:,:,r)); caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T1NoMFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(3); clf;
imagesc(T1_map_masked_MFI(:,:,r)); caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T1MFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(4); clf;
imagesc(T2_map_masked(:,:,r)); caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T2NoMFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(5); clf;
imagesc(T2_map_masked_MFI(:,:,r)); caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T2MFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

% row 3
r = 4;

figure(1); clf;
imagesc(B0_map_masked(:,:,r));
caxis([-450 450])
set(gcf,'Color','w');
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(25,190,'c','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('imgB0_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc(T1_map_masked(:,:,r)); caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T1NoMFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(3); clf;
imagesc(T1_map_masked_MFI(:,:,r)); caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T1MFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(4); clf;
imagesc(T2_map_masked(:,:,r)); caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T2NoMFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(5); clf;
imagesc(T2_map_masked_MFI(:,:,r)); caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T2MFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

% row 4
r = 3;

figure(1); clf;
imagesc(B0_map_masked(:,:,r));
caxis([-450 450])
set(gcf,'Color','w');
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(25,190,'d','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('imgB0_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc(T1_map_masked(:,:,r)); caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
c = colorbar('south','FontSize',fontSzCBTicks,'LineWidth',2.0);
c.Label.String = 'T1 (ms)';
c.Label.FontSize = fontSzCBlab;
c.Color = 'w';
c.Position(1) = 0.30;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T1NoMFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(3); clf;
imagesc(T1_map_masked_MFI(:,:,r)); caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
c = colorbar('south','FontSize',fontSzCBTicks,'LineWidth',2.0);
c.Label.String = 'T1 (ms)';
c.Label.FontSize = fontSzCBlab;
c.Color = 'w';
c.Position(1) = 0.30;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T1MFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(4); clf;
imagesc(T2_map_masked(:,:,r)); caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
c = colorbar('south','FontSize',fontSzCBTicks,'LineWidth',2.0);
c.Label.String = 'T2 (ms)';
c.Label.FontSize = fontSzCBlab;
c.Color = 'w';
c.Position(1) = 0.30;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T2NoMFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

figure(5); clf;
imagesc(T2_map_masked_MFI(:,:,r)); caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
c = colorbar('south','FontSize',fontSzCBTicks,'LineWidth',2.0);
c.Label.String = 'T2 (ms)';
c.Label.FontSize = fontSzCBlab;
c.Color = 'w';
c.Position(1) = 0.30;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
txt = sprintf('T2MFI_%d = imread([fig_out ''tmp.png'']);',r);
eval(txt);
delete([fig_out 'tmp.png'])

% form total figure
r_v = [1,2,4,3];
nrows = 4;
fig2_m = [];
for ii = 1:nrows
    txt = sprintf('fig2_m = cat(1,fig2_m,[repmat(imgB0_%d,[1 1 3]),T1NoMFI_%d,T1MFI_%d,T2NoMFI_%d,T2MFI_%d]);',r_v(ii)*ones(1,5));
    eval(txt);
end

% add divider between B0 and T1 and T2 and between rows
border_width1 = 32;
border_width2 = 64;
[rtmp,ctmp,~] = size(fig2_m);
[rmap,cmap,~] = size(imgB0_1);
fig2_m(:,(cmap + 1):(cmap + border_width1),:) = ...
    255*ones(rtmp,border_width1,3);
fig2_m(:,((3*cmap)-border_width1/2+1):((3*cmap)+border_width1/2),:) = ...
    255*ones(rtmp,border_width1,3);
fig2_m(rmap-border_width2/2 + 1:rmap+border_width2/2,:,:) = ...
    255*ones(border_width2,ctmp,3);
fig2_m(2*rmap-border_width2/2 + 1:2*rmap+border_width2/2,:,:) = ...
    255*ones(border_width2,ctmp,3);
fig2_m(3*rmap-border_width2/2 + 1:3*rmap+border_width2/2,:,:) = ...
    255*ones(border_width2,ctmp,3);

figure(10); clf;
imshow(fig2_m)
txtexec = sprintf('export_fig([fig_out ''MRI_Figure2_MRFMFI.png''],''-r%d'')',resImg);
eval(txtexec);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure2_MRFMFI.tif''],''-r%d'')',resImg);
eval(txtexec);

%% Figure 3 MRI system phantom load data 

close all

% files to be loaded
fn1 = 'MRF_proc_HPD_top_B0het.mat';
fn1MFI = 'MRF_proc_HPD_top_B0het_MFI.mat';
fn1B0 = 'MRF_B0_HPD_top_B0het.mat';

fn2 = 'MRF_proc_HPD_mid_B0het.mat';
fn2MFI = 'MRF_proc_HPD_mid_B0het_MFI.mat';
fn2B0 = 'MRF_B0_HPD_mid_B0het.mat';

img_mag = zeros(effMtx,effMtx,2);
T1_map_masked = img_mag;
T2_map_masked = img_mag;
img_mag_MFI = img_mag;
T1_map_masked_MFI = img_mag;
T2_map_masked_MFI = img_mag;
B0_map_masked = img_mag;

for ii = 1:2

    txt = sprintf('load([datadir fn%d])',ii);
    eval(txt);
    
    my_img = abs(sum(input_MRF_dp.img_AR_comb,3));
    img_mag(:,:,ii) = my_img;
    mythresh_v = multithresh(my_img,MFI_thresh);
    my_T1_map = output_MRF_match.T1_map;
    my_T2_map = output_MRF_match.T2_map;
    my_T1_map(my_img < mythresh_v(1)) = 0;
    my_T2_map(my_img < mythresh_v(1)) = 0;
    T1_map_masked(:,:,ii) = my_T1_map;
    T2_map_masked(:,:,ii) = my_T2_map;
    
    txt = sprintf('load([datadir fn%dMFI])',ii);
    eval(txt);
    txt = sprintf('load([datadir_B0 fn%dB0])',ii);
    eval(txt);
    
    my_img = abs(sum(input_MRF_dp.img_AR_comb,3));
    img_mag_MFI(:,:,ii) = my_img;
    mythresh_v = multithresh(my_img,MFI_thresh);
    my_T1_map = output_MRF_match.T1_map;
    my_T2_map = output_MRF_match.T2_map;
    my_T1_map(my_img < mythresh_v(1)) = 0;
    my_T2_map(my_img < mythresh_v(1)) = 0;
    T1_map_masked_MFI(:,:,ii) = my_T1_map;
    T2_map_masked_MFI(:,:,ii) = my_T2_map;
    B0_map_masked(:,:,ii) = output_B0.B0_map_masked;

end

%% Figure 3 MRI system phantom create figure

% desired image resolution in dpi
resImg = 500;

% font sizes
fontSzFigLet = 40; % figure letter size
fontSzCBlab = 28; % colorbar label
fontSzCBTicks = 28; % colorbar tick labels

figure(1); clf;
imagesc(B0_map_masked(:,:,1));
set(gcf,'Color','w');
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
c = colorbar('FontSize',fontSzCBTicks,'LineWidth',2.0);
c.Label.String = 'B0 off-resonance (Hz)';
c.Label.FontSize = fontSzCBlab;
c.Position(1) = 0.83;
c.Position(2) = 0.15;
c.Position(3) = 0.03;
c.Position(4) = 0.75;
text(25,240-25,'c','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')

txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgB0_1 = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc(img_mag(:,:,1));
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(50,25,'w/o MFI','FontSize',2*fontSzFigLet,'Color','w','FontWeight','bold')
text(25,240-25,'a','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgNoMFI_1 = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(3); clf;
imagesc(img_mag_MFI(:,:,1));
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(85,25,'MFI','FontSize',2*fontSzFigLet,'Color','w','FontWeight','bold')
text(25,240-25,'b','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgMFI_1 = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(1); clf;
imagesc(B0_map_masked(:,:,2));
set(gcf,'Color','w');
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
c = colorbar('FontSize',fontSzCBTicks,'LineWidth',2.0);
c.Label.String = 'B0 off-resonance (Hz)';
c.Label.FontSize = fontSzCBlab;
c.Position(1) = 0.83;
c.Position(2) = 0.15;
c.Position(3) = 0.03;
c.Position(4) = 0.75;
text(25,240-25,'f','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgB0_2 = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc(img_mag(:,:,2));
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(25,240-25,'d','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgNoMFI_2 = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(3); clf;
imagesc(img_mag_MFI(:,:,2));
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(25,240-25,'e','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgMFI_2 = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

imgFig3 = [imgNoMFI_1 imgMFI_1 imgB0_1; imgNoMFI_2 imgMFI_2 imgB0_2];
figure(4); clf;
imshow(imgFig3)
txtexec = sprintf('export_fig([fig_out ''MRI_Figure3_MRFMFI.png''],''-r%d'')',resImg);
eval(txtexec);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure3_MRFMFI.tif''],''-r%d'')',resImg);
eval(txtexec);

%% Figure 4 MRI system T1 and T2 maps

close all

% desired image resolution in dpi
resImg = 500;

% font sizes
fontSzFigLet = 40; % figure letter size
fontSzCBlab = 28; % colorbar label
fontSzCBTicks = 28; % colorbar tick labels

minT1 = 0;
maxT1 = 2000;
minT2 = 0;
maxT2 = 550;

figure(1); clf;
imagesc(T1_map_masked(:,:,1));
caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
text(50,25,'w/o MFI','FontSize',2*fontSzFigLet,'Color','w','FontWeight','bold')
text(25,240-25,'a','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
T1_noMFI_m = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc([T1_map_masked_MFI(:,:,1) zeros(effMtx,50)]);
caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
text(85,25,'MFI','FontSize',2*fontSzFigLet,'Color','w','FontWeight','bold')
text(25,240-25,'b','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
c = colorbar('FontSize',fontSzCBTicks,'LineWidth',2.0,'Color','w');
c.Label.String = 'T1 (ms)';
c.Label.FontSize = fontSzCBlab;
c.Position(1) = 0.73;
c.Position(2) = 0.15;
c.Position(3) = 0.03;
c.Position(4) = 0.75;
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
T1_MFI_m = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

% T2

figure(1); clf;
imagesc(T2_map_masked(:,:,2));
caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
text(25,240-25,'c','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
T2_noMFI_m = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc([T2_map_masked_MFI(:,:,2) zeros(effMtx,50)]);
caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
text(25,240-25,'d','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
c = colorbar('FontSize',fontSzCBTicks,'LineWidth',2.0,'Color','w');
c.Label.String = 'T2 (ms)';
c.Label.FontSize = fontSzCBlab;
c.Position(1) = 0.73;
c.Position(2) = 0.15;
c.Position(3) = 0.03;
c.Position(4) = 0.75;
eval(txtexec)
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
T2_MFI_m = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

fig4_m = [T1_noMFI_m T1_MFI_m; T2_noMFI_m T2_MFI_m];

figure(3); clf;
imshow(fig4_m);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure4_MRFMFI.png''],''-r%d'')',resImg);
eval(txtexec);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure4_MRFMFI.tif''],''-r%d'')',resImg);
eval(txtexec);

%% Table 1 MRI system quant results (T1)

close all

% filename
fn1_stats = 'ROIstats_HPD_2017.mat';

% get HPD statistics
analyze_HPD

nroi = 14;
maxvox = 5000; % maximum possible number of voxels in an roi

% load data
load([datadir fn1_stats])

% T1 boxplots
X_m = zeros(maxvox,nroi*4);
HPD_v = zeros(1,nroi*4);
G1_c = zeros(1,nroi*4);
G2_c = cell(1,nroi*4);
n = 1;

for t = 1:1
    for kk = 1:nroi
        for mfi = 1:4
            
            % get measurements
            txt = sprintf('x_v = roistats.d%d.sl%d.T%d.msr%d;',mfi,t,t,kk);
            eval(txt);
            nx = numel(x_v);
            X_m(1:nx,n) = x_v;
            X_m(nx+1:end,n) = NaN;
            txt = sprintf('tmp_v = roistats.HPD.sl%d.T%d.mean;',t,t);
            eval(txt);
            HPD_v(n) = tmp_v(kk);
            
            % make labels
            G1_c(n) = kk;
            switch mfi
                case 1
                    g2 = 'ws'; %'U no MFI';
                case 2
                    g2 = 'ws+'; %'U MFI';
                case 3
                    g2 = 'ps'; %'NU no MFI';
                case 4
                    g2 = 'ps+'; %'NU MFI';
            end
            
            G2_c(n) = {g2};
            
            n = n + 1;
        end
    end
end

    
% compute stats
% coefficient of variation
X_T1stats_m = zeros(4,n-1);
for ii = 1:n-1
    tmp_v = X_m(~isnan(X_m(:,ii)),ii);
    X_T1stats_m(1,ii) = std(tmp_v);
    X_T1stats_m(2,ii) = mean(tmp_v);
    X_T1stats_m(4,ii) = median(tmp_v);
end

X_T1stats_m(3,:) = X_T1stats_m(1,:)./X_T1stats_m(2,:); % coeff of variation

T1_coeffvar_m = (reshape(X_T1stats_m(3,:),[4 (n-1)/4]))';

% CCC with CIs
T1_CCC_m = zeros(4,3);
for ii = 1:4
    idx_v = ii:4:(n-1);
    tmp_msr_v = X_T1stats_m(2,idx_v);
    tmp_HPD_v = HPD_v(idx_v);
    [pc,ci_top,cb,u,pr] = ccc(tmp_HPD_v,tmp_msr_v);
    T1_CCC_m(ii,:) = [pc ci_top];
end

T_T1 = array2table(T1_CCC_m,...
    'RowNames',{'T1_ws','T1_ws+','T1_ps','T1_ps+'},...
    'VariableNames',{'Mean','Lower_95_CI','Upper_95_CI'})

%% Figure 5 MRI system T1 boxplots

figure(1); clf;
set(gcf,'Color','w','Position',[100 200 800 600]);
boxplot(X_m(:,1:7*4),{G1_c(1:7*4), G2_c(1:7*4)},'colors',repmat('rgbm',1,nroi/2),'plotstyle','compact')
ylabel('T1 (ms)')
set(gca,'FontSize',fontSzCBlab,'LineWidth',2)
h = findobj(gca,'Type','text');
for ii = 1:4*nroi
    tmp = h(ii);
    tmp.FontSize = fontSzCBTicks;
    if ii <= 4*nroi/2
        tmp.Position(2) = tmp.Position(2) - 20;
        tmp.Rotation = 90;
    else
        tmp.Position(2) = tmp.Position(2) - 25;
        tmp.Position(1) = tmp.Position(1) + 20;
        tmp.Rotation = 0;
    end
end

% add HPD reference lines
xstep = 4.6;
idx_v = 1:4:4*nroi/2;
for ii = 1:nroi/2
    a = line([(0 + xstep*(ii-1)) ii*xstep],[HPD_v(idx_v(ii)) HPD_v(idx_v(ii))]);
    a.LineStyle = '--';
    a.LineWidth = 3;
    a.Color = 'k';
end

text(2,250,'a','FontSize',fontSzFigLet,'Color','k','FontWeight','bold');

txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
Fig5a = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

% 5b

figure(2); clf;
set(gcf,'Color','w','Position',[100 200 800 600]);
boxplot(X_m(:,((7*4)+1):end),{G1_c(7*4+1:end), G2_c(7*4+1:end)},'colors',repmat('rgbm',1,nroi/2),'plotstyle','compact')
ylabel('T1 (ms)')
set(gca,'FontSize',fontSzCBlab,'LineWidth',2)
h = findobj(gca,'Type','text');
for ii = 1:4*nroi
    tmp = h(ii);
    tmp.FontSize = fontSzCBTicks;
    if ii <= 4*nroi/2
        tmp.Position(2) = tmp.Position(2) - 20;
        tmp.Rotation = 90;
    else
        tmp.Position(2) = tmp.Position(2) - 25;
        tmp.Position(1) = tmp.Position(1) + 20;
        tmp.Rotation = 0;
    end
end

% add HPD reference lines
xstep = 4.6;
idx_v = 4*nroi/2 + 1:4:4*nroi;
for ii = 1:nroi/2
    a = line([(0 + xstep*(ii-1)) ii*xstep],[HPD_v(idx_v(ii)) HPD_v(idx_v(ii))]);
    a.LineStyle = '--';
    a.LineWidth = 3;
    a.Color = 'k';
end

text(2,35,'b','FontSize',fontSzFigLet,'Color','k','FontWeight','bold');

txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
Fig5b = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

%% Table 1 MRI system quant results (T2)

X_m = zeros(maxvox,nroi*4);
HPD_v = zeros(1,nroi*4);
G1_c = zeros(1,nroi*4);
G2_c = cell(1,nroi*4);
n = 1;

for t = 2:2
    for kk = 1:nroi
        for mfi = 1:4
            
            % get measurements
            txt = sprintf('x_v = roistats.d%d.sl%d.T%d.msr%d;',mfi,t,t,kk);
            eval(txt);
            nx = numel(x_v);
            X_m(1:nx,n) = x_v;
            X_m(nx+1:end,n) = NaN;
            txt = sprintf('tmp_v = roistats.HPD.sl%d.T%d.mean;',t,t);
            eval(txt);
            HPD_v(n) = tmp_v(kk);
            
            % make labels
            G1_c(n) = kk;
            switch mfi
                case 1
                    g2 = 'ws'; %'U no MFI';
                case 2
                    g2 = 'ws+'; %'U MFI';
                case 3
                    g2 = 'ps'; %'NU no MFI';
                case 4
                    g2 = 'ps+'; %'NU MFI';
            end
            
            G2_c(n) = {g2};
            
            n = n + 1;
        end
    end
end

% compute stats

X_T2stats_m = zeros(4,n-1);
for ii = 1:n-1
    tmp_v = X_m(~isnan(X_m(:,ii)),ii);
    X_T2stats_m(1,ii) = std(tmp_v);
    X_T2stats_m(2,ii) = mean(tmp_v);
    X_T2stats_m(4,ii) = median(tmp_v);
end

X_T2stats_m(3,:) = X_T2stats_m(1,:)./X_T2stats_m(2,:); % coeff of variation

T2_coeffvar_m = (reshape(X_T2stats_m(3,:),[4 (n-1)/4]))';

% CCC with CIs
T2_CCC_m = zeros(4,3);
for ii = 1:4
    idx_v = ii:4:(n-1);
    tmp_msr_v = X_T2stats_m(2,idx_v);
    tmp_HPD_v = HPD_v(idx_v);
    [pc,ci_top,cb,u,pr] = ccc(tmp_HPD_v,tmp_msr_v);
    T2_CCC_m(ii,:) = [pc ci_top];
end

T_T2 = array2table(T2_CCC_m,...
    'RowNames',{'T2_ws','T2_ws+','T2_ps','T2_ps+'},...
    'VariableNames',{'Mean','Lower_95_CI','Upper_95_CI'})

%% Figure 5 MRI system T2 boxplots

figure(3); clf;
set(gcf,'Color','w','Position',[100 200 800 600]);
boxplot(X_m(:,1:7*4),{G1_c(1:7*4), G2_c(1:7*4)},'colors',repmat('rgbm',1,nroi/2),'plotstyle','compact')
ylabel('T2 (ms)')
set(gca,'FontSize',fontSzCBlab,'LineWidth',2)
h = findobj(gca,'Type','text');
for ii = 1:4*nroi
    tmp = h(ii);
    tmp.FontSize = fontSzCBTicks;
    if ii <= 4*nroi/2
        tmp.Position(2) = tmp.Position(2) - 20;
        tmp.Rotation = 90;
    else
        tmp.Position(2) = tmp.Position(2) - 25;
        tmp.Position(1) = tmp.Position(1) + 20;
        tmp.Rotation = 0;
    end
end

% add HPD reference lines
xstep = 4.6;
idx_v = 1:4:4*nroi/2;
for ii = 1:nroi/2
    a = line([(0 + xstep*(ii-1)) ii*xstep],[HPD_v(idx_v(ii)) HPD_v(idx_v(ii))]);
    a.LineStyle = '--';
    a.LineWidth = 3;
    a.Color = 'k';
end

text(2,75,'c','FontSize',fontSzFigLet,'Color','k','FontWeight','bold');

txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
Fig5c = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

% 5d

figure(4); clf;
set(gcf,'Color','w','Position',[100 200 800 600]);
boxplot(X_m(:,((7*4)+1):end),{G1_c(7*4+1:end), G2_c(7*4+1:end)},'colors',repmat('rgbm',1,nroi/2),'plotstyle','compact');
ylabel('T2 (ms)')
set(gca,'FontSize',fontSzCBlab,'LineWidth',2)
h = findobj(gca,'Type','text');
for ii = 1:4*nroi
    tmp = h(ii);
    tmp.FontSize = fontSzCBTicks;
    if ii <= 4*nroi/2
        tmp.Position(2) = tmp.Position(2) - 20;
        tmp.Rotation = 90;
    else
        tmp.Position(2) = tmp.Position(2) - 25;
        tmp.Position(1) = tmp.Position(1) + 20;
        tmp.Rotation = 0;
    end
end

% add HPD reference lines
xstep = 4.6;
idx_v = 4*nroi/2 + 1:4:4*nroi;
for ii = 1:nroi/2
    a = line([(0 + xstep*(ii-1)) ii*xstep],[HPD_v(idx_v(ii)) HPD_v(idx_v(ii))]);
    a.LineStyle = '--';
    a.LineWidth = 3;
    a.Color = 'k';
end

text(2,18,'d','FontSize',fontSzFigLet,'Color','k','FontWeight','bold');

txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
Fig5d = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

[r4a,c4a,~] = size(Fig5a);
[r4b,c4b,~] = size(Fig5b);
[r4c,c4c,~] = size(Fig5c);
[r4d,c4d,~] = size(Fig5d);

cxtra = c4a-c4c;

imgFig5 = [[Fig5a; 255*ones(1,c4a,3)] Fig5b; ...
    [Fig5c; 255*ones(1,c4c,3)] 255*ones(r4d,cxtra,3) Fig5d];
figure(5); clf;
imshow(imgFig5)
txtexec = sprintf('export_fig([fig_out ''MRI_Figure5_MRFMFI.png''],''-r%d'')',resImg);
eval(txtexec);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure5_MRFMFI.tif''],''-r%d'')',resImg);
eval(txtexec);

%% Figure 6 coefficients of variation

close all

% desired image resolution in dpi
resImg = 500;

% font sizes
fontSzFigLet = 40; % figure letter size
fontSzCBlab = 28; % colorbar label
fontSzCBTicks = 28; % colorbar tick labels

figure(1); clf;
barh = bar(T1_coeffvar_m*100,'LineWidth',2);
barh(1).FaceColor = 'r' ;
barh(2).FaceColor = 'g';
barh(3).FaceColor = 'b';
barh(4).FaceColor = 'm';
xlabel('ROI')
ylabel('Coefficient of variation (%)')
set(gcf,'Color','w','Position',[100 200 800 600]);
set(gca,'LineWidth',2,'FontSize',fontSzCBlab)
l = legend('ws','ws+','ps','ps+');
l.Position(1) = 0.65;
l.Position(2) = 0.70;
text(2,100,'T1','FontSize',fontSzFigLet*2,'Color','k','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
T1cov_m = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(2); clf;
barh = bar(T2_coeffvar_m*100,'LineWidth',2);
barh(1).FaceColor = 'r';
barh(2).FaceColor = 'g';
barh(3).FaceColor = 'b';
barh(4).FaceColor = 'm';
xlabel('ROI')
ylabel('Coefficient of variation (%)')
set(gcf,'Color','w','Position',[100 200 800 600]);
set(gca,'LineWidth',2,'FontSize',fontSzCBlab)
l = legend('ws','ws+','ps','ps+');
l.Position(1) = 0.65;
l.Position(2) = 0.70;
text(2,83,'T2','FontSize',fontSzFigLet*2,'Color','k','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
T2cov_m = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

%% handle possible slight mismatch of PNG sizes
png_size_1 = min( size(T1cov_m,1), size(T2cov_m,1));
imgFig6 = [T1cov_m(1:png_size_1,:,:) 255*ones(png_size_1,20,3) T2cov_m(1:png_size_1,:,:) ];

figure(3); clf;
imshow(imgFig6)
txtexec = sprintf('export_fig([fig_out ''MRI_Figure6_MRFMFI.png''],''-r%d'')',resImg);
eval(txtexec);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure6_MRFMFI.tif''],''-r%d'')',resImg);
eval(txtexec);

% get legend and add to Figure 5
mrsys_legend_m = imgFig6(250:1040,8415:9055,:);
[nrleg,ncleg,~] = size(mrsys_legend_m);
imgFig5(770:770+nrleg-1,8200:8200+ncleg-1,:) = mrsys_legend_m;
imgFig5(770+3650:770+3650+nrleg-1,8200:8200+ncleg-1,:) = mrsys_legend_m;
figure(5); clf;
imshow(imgFig5)
txtexec = sprintf('export_fig([fig_out ''MRI_Figure5_MRFMFI.png''],''-r%d'')',resImg);
eval(txtexec);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure5_MRFMFI.tif''],''-r%d'')',resImg);
eval(txtexec);


%% Figure 7 in vivo brain load data

close all

% files to be loaded
fn1 = 'MRF_proc_invivo_sl2_long.mat';
fn1MFI = 'MRF_proc_invivo_sl2_long_MFI.mat';
fn1B0 = 'MRF_B0_invivo_sl2.mat';
fn1T1W = 'T1W_invivo.mat';
fn1mask = 'T1W_invivo_brain_mask.mat';

load([datadir_B0 fn1mask]);
load([datadir fn1]);
load([datadir_B0 fn1B0]);
load([datadir_B0 fn1T1W]);

my_thresh = multithresh(img_T1W,1);

T1_map_masked = output_MRF_match.T1_map;
T1_map_masked(~brain_mask) = 0;
T2_map_masked = output_MRF_match.T2_map;
T2_map_masked(~brain_mask) = 0;

load([datadir fn1MFI]);
T1_map_masked_MFI = output_MRF_match.T1_map;
T1_map_masked_MFI(~brain_mask) = 0;
T2_map_masked_MFI = output_MRF_match.T2_map;
T2_map_masked_MFI(~brain_mask) = 0;
B0_map_masked = output_B0.B0_map_masked;

% for signal plots
ri_CSF = 63; ci_CSF = 124;
ri_WM = 160; ci_WM = 158;

%% Figure 7 in vivo brain

% desired image resolution in dpi
resImg = 500;

% font sizes
fontSzFigLet = 40; % figure letter size
fontSzCBlab = 28; % colorbar label
fontSzCBTicks = 24; % colorbar tick labels

minT1 = 0;
maxT1 = 2500;
minT2 = 0;
maxT2 = 300;

figure(1); clf;
imagesc(B0_map_masked);
set(gcf,'Color','w');
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
caxis([-150 150])
c = colorbar('FontSize',fontSzCBTicks,'LineWidth',2.0,'Color','w');
c.Label.String = 'B0 off-resonance (Hz)';
c.Label.FontSize = fontSzCBlab-2;
c.Label.Position = [-3.2 0 0];
c.Position(1) = 0.31;
c.Position(2) = 0.15;
c.Position(3) = 0.02;
c.Position(4) = 0.75;
text(55,25,'a','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgB0 = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(2); clf;
imagesc(img_T1W);
set(gcf,'Color','w');
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(gray);
text(25,25,'d','FontSize',fontSzFigLet,'Color','w','FontWeight','bold')
hold on
plot(ci_CSF,ri_CSF,'wx','LineWidth',6,'MarkerSize',20)
plot(ci_WM,ri_WM,'w+','LineWidth',6,'MarkerSize',20)
hold off
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
imgT1w = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(3); clf;
imagesc([T1_map_masked; T1_map_masked_MFI]); caxis([minT1 maxT1])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
c = colorbar('FontSize',fontSzCBTicks/2,'LineWidth',2.0);
c.Color = 'w';
c.Label.String = 'T1 (ms)';
c.Label.FontSize = fontSzCBlab/2;
c.Label.Position = [-1.8 (maxT1-minT1)/2 0];
c.Position(1) = 0.62;
c.Position(2) = 0.3;
c.Position(3) = 0.01;
c.Position(4) = 0.4;
text(25,25,'b','FontSize',fontSzFigLet/2,'Color','w','FontWeight','bold')
text(25,265,'e','FontSize',fontSzFigLet/2,'Color','w','FontWeight','bold')
annotation('arrow',[0.51 0.51],[0.805 0.805],'Color','w','HeadLength',10,'Headwidth',10)
annotation('arrow',[0.51 0.51],[0.395 0.395],'Color','w','HeadLength',10,'Headwidth',10)
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',2*resImg);
eval(txtexec)
T1map = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

figure(4); clf;
imagesc([T2_map_masked; T2_map_masked_MFI]); caxis([minT2 maxT2])
axis image
axis off
set(gcf,'Color','w','Position',[100 200 800 600]);
colormap(IDLrainbow2);
c = colorbar('FontSize',fontSzCBTicks/2,'LineWidth',2.0);
c.Color = 'w';
c.Label.String = 'T2 (ms)';
c.Label.FontSize = fontSzCBlab/2;
c.Label.Position = [-1.8 (maxT2-minT2)/2 0];
c.Position(1) = 0.62;
c.Position(2) = 0.3;
c.Position(3) = 0.01;
c.Position(4) = 0.4;
text(25,25,'c','FontSize',fontSzFigLet/2,'Color','w','FontWeight','bold')
text(25,265,'f','FontSize',fontSzFigLet/2,'Color','w','FontWeight','bold')
annotation('arrow',[0.51 0.51],[0.82 0.82],'Color','w','HeadLength',10,'Headwidth',10)
annotation('arrow',[0.51 0.51],[0.41 0.41],'Color','w','HeadLength',10,'Headwidth',10)
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',2*resImg);
eval(txtexec)
T2map = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

[nrB0,ncB0,~] = size(imgB0);
[nrT1map,ncT1map,~] = size(T1map);

cxtra = 2*ncB0 - 2*ncT1map;

tmp_m = repmat([imgB0; imgT1w],[1 1 3]);
rxtra = size(tmp_m,1) - nrT1map;

imgFig7 = [tmp_m [T1map T2map; zeros(rxtra,2*ncT1map,3)]]; % horizontal

figure(5); clf;
imshow(imgFig7)
txtexec = sprintf('export_fig([fig_out ''MRI_Figure7_MRFMFI.png''],''-r%d'')',resImg);
eval(txtexec);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure7_MRFMFI.tif''],''-r%d'')',resImg);
eval(txtexec);


%% Fig 8 signal evolutions from in vivo brain

close all

% desired image resolution in dpi
resImg = 500;

% font sizes
fontSzFigLet = 40; % figure letter size
fontSzCBlab = 28; % colorbar label
fontSzCBTicks = 24; % colorbar tick labels

lw = 2; % line width

load([datadir fn1]);
s_v = squeeze(input_MRF_dp.img_AR_comb(ri_CSF,ci_CSF,:));
load([datadir fn1MFI]);
s_MFI_v = squeeze(input_MRF_dp.img_AR_comb(ri_CSF,ci_CSF,:));

s_v = s_v./norm(s_v);
s_MFI_v = s_MFI_v./norm(s_MFI_v);

figure(1); clf;
set(gcf,'Color','w','Position',[100 200 800 600]);
plot((abs(s_v)),'LineWidth',lw)
hold on
plot((abs(s_MFI_v)),'--','LineWidth',lw)
hold off
set(gca,'fontsize',fontSzCBTicks)
ylim([0 0.10])
T1 = T1_map_masked(ri_CSF,ci_CSF);
T2 = T2_map_masked(ri_CSF,ci_CSF);
T1_MFI = T1_map_masked_MFI(ri_CSF,ci_CSF);
T2_MFI = T2_map_masked_MFI(ri_CSF,ci_CSF);
B0_Hz = B0_map_masked(ri_CSF,ci_CSF);
h = legend(['w/o MFI - T1 ' num2str(T1) ' ms - T2 ' num2str(T2) ' ms'],...
    ['MFI - T1 ' num2str(T1_MFI) ' ms - T2 ' num2str(T2_MFI) ' ms']);
h.Position(1) = 0.27;
h.Position(2) = 0.84;
xlabel('TR')
ylabel('normalized signal')
text(300,0.082,['B0 = ' num2str(round(B0_Hz)) ' Hz'],'FontSize',2*fontSzCBTicks,'Color','k','FontWeight','bold')
text(25,0.09,'a','FontSize',fontSzFigLet,'Color','k','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
CSF_plot = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

load([datadir fn1]);
s_v = squeeze(input_MRF_dp.img_AR_comb(ri_WM,ci_WM,:));
load([datadir fn1MFI]);
s_MFI_v = squeeze(input_MRF_dp.img_AR_comb(ri_WM,ci_WM,:));

s_v = s_v./norm(s_v);
s_MFI_v = s_MFI_v./norm(s_MFI_v);

figure(2); clf;
set(gcf,'Color','w','Position',[100 200 800 600]);
plot((abs(s_v)),'LineWidth',lw)
hold on
plot((abs(s_MFI_v)),'--','LineWidth',lw)
hold off
set(gca,'fontsize',fontSzCBTicks)
ylim([0 0.10])
T1 = T1_map_masked(ri_WM,ci_WM);
T2 = T2_map_masked(ri_WM,ci_WM);
T1_MFI = T1_map_masked_MFI(ri_WM,ci_WM);
T2_MFI = T2_map_masked_MFI(ri_WM,ci_WM);
B0_Hz = B0_map_masked(ri_WM,ci_WM);
h = legend(['w/o MFI - T1 ' num2str(T1) ' ms - T2 ' num2str(T2) ' ms'],...
    ['MFI - T1 ' num2str(T1_MFI) ' ms - T2 ' num2str(T2_MFI) ' ms']);
h.Position(1) = 0.27;
h.Position(2) = 0.84;
xlabel('TR')
ylabel('normalized signal')
text(300,0.082,['B0 = ' num2str(round(B0_Hz)) ' Hz'],'FontSize',2*fontSzCBTicks,'Color','k','FontWeight','bold')
text(25,0.09,'b','FontSize',fontSzFigLet,'Color','k','FontWeight','bold')
txtexec = sprintf('export_fig([fig_out ''tmp.png''],''-r%d'')',resImg);
eval(txtexec)
WM_plot = imread([fig_out 'tmp.png']);
delete([fig_out 'tmp.png'])

imgFig8 = [CSF_plot; WM_plot];

figure(5); clf;
imshow(imgFig8)
txtexec = sprintf('export_fig([fig_out ''MRI_Figure8_MRFMFI.png''],''-r%d'')',resImg);
eval(txtexec);
txtexec = sprintf('export_fig([fig_out ''MRI_Figure8_MRFMFI.tif''],''-r%d'')',resImg);
eval(txtexec);
