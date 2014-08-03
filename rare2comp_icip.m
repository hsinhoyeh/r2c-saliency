function [ salmap ] = ...
    rare2comp_icip( filenames)
%idea: find salient sopts and spread it to compactness
%spatial wrighting

ratio = 0.05;
verbose=0;
%0 for viewing nothing
%1 for displaying the processing message
%2 for view the salient result only


salmap = cell(length(filenames), 1);
for ii=1: length(filenames)
    ori_im = im2double(imread( filenames{ii,1} ));
    [ori_h, ori_w, ori_c] = size(ori_im);
    
    %Lab color space
    [im1, im2, im3]=RGB2Lab(ori_im);
    %additional matlab file of RGB2Lab
    
    fim(:,:,1)=mat2gray(im1);%feature image
    fim(:,:,2)=mat2gray(im2);
    fim(:,:,3)=mat2gray(im3);

    %%
    %generate finer layer feature representation
    fparams=[];
    fparams.pwidth=2;
    [fpmat, fpcmat, foparams] = ...%finer patches
        patch_gen(fim, fparams);

    fnum_patch=foparams.num_patchx*foparams.num_patchy;
    fpmatim = reshape([1:1: fnum_patch], ...
        foparams.num_patchy, foparams.num_patchx);
    %pmatim is the feature image

    
    %absolute contrast
    mean_vec = mean(fpmat,1);
    acdifmat= fpmat - repmat(mean_vec, fnum_patch,1);
    acmat = sqrt(sum(acdifmat.^2,2));
    acdifmat=[];
    %bin contrast

    tri_fp2fpdist=pdist(fpmat, 'euclidean');
    squ_fp2fpdist=squareform(tri_fp2fpdist);
    num_bin = 2; %%tunning parameter
    bin_sigma=1;  %%tunning parameter
    [ncutlabel, ~, ~] = ...
        ncutW(exp(-squ_fp2fpdist), num_bin);
    binlabel=ncutlabel*[1:1:num_bin]';
    
    bincnt=histc(binlabel, [1:1:num_bin]);
    binrate=bincnt ./ sum(bincnt);
    bcexp = exp(-binrate/ bin_sigma );
    bcmat=bcexp(binlabel);
    
    binlabel=[];
    %local contrast
    lc_lambda=1;
    tri_spa2spadist=pdist( fpcmat(:,[1,2]), 'euclidean' );
    sqt_spa2spadist=squareform(tri_spa2spadist);
    tri_lcdist=tri_fp2fpdist ./ (1+lc_lambda*tri_spa2spadist);
    lcmat = mean(squareform(tri_lcdist),2);
    
    tri_lcdist=[];
    
    gspa=[0.5,0.5];
    gsloc=[ori_w*gspa(:,1), ori_h*gspa(:,2)];
    
    cgdist=pdist2(fpcmat(:,[1:2]),gsloc, 'euclidean');
    %coarse patch to golden patch distmat
    cmindist=min(cgdist,[],2);%coarse minimum distance
    spasigma=80;%for spatial weighting
    cmprob=exp(-cmindist/spasigma);

    
    fusemat = cmprob.*acmat.*bcmat.*lcmat;

    %%
    %discover sps
    spot_rate=ratio;
    spot_rth=ceil(spot_rate*fnum_patch);
    sortfumat=sort(fusemat, 'descend');
    spot_th=sortfumat(spot_rth,1);
    
    spindmat=fusemat >= spot_th;
    spimat=spindmat.*fusemat;
    %spread it
    sp_delth=0.01;%for spatial
    sp_sigma=0.5;%for color
    
    sp_colmat = squ_fp2fpdist(spindmat,:);
    %salient spots for color projection
    sp_spamat = sqt_spa2spadist(spindmat,:);
    %salient spots for spatial projection
    sp_salmat = fusemat(spindmat,1);
    sp_salmat = repmat(sp_salmat, 1, fnum_patch);
    %salient spots for saliency projection
    t1= exp(-sp_sigma*sp_colmat);%color compactness
    t2= exp(-sp_delth*sp_spamat);%spatial compactness 
    c1 = max(t1,[],1)';
    c2 = max(t2,[],1)';
    compness=c1.*c2;
    
    spim = zeros(ori_h, ori_w);
    for pp=1:fnum_patch
        aa=fpcmat(pp,[3,4]);
        bb=fpcmat(pp,[5,6]);
        spim(aa(2):bb(2), aa(1):bb(1)) = ...
            compness(pp,1); 
    end

if verbose >= 2

    %%for visualization
    acim = zeros(ori_h, ori_w);
    bcim = zeros(ori_h, ori_w);
    lcim = zeros(ori_h, ori_w);
    gwim = zeros(ori_h, ori_w);
    c1im = zeros(ori_h, ori_w);
    c2im = zeros(ori_h, ori_w);
    spiim = zeros(ori_h, ori_w);
    fuim = zeros(ori_h, ori_w);
    
    for pp=1:fnum_patch
        aa=fpcmat(pp,[3,4]);
        bb=fpcmat(pp,[5,6]);
        acim(aa(2):bb(2), aa(1):bb(1)) = ...
            acmat(pp,1);
        bcim(aa(2):bb(2), aa(1):bb(1)) = ...
            bcmat(pp,1);
        lcim(aa(2):bb(2), aa(1):bb(1)) = ...
            lcmat(pp,1);
        gwim(aa(2):bb(2), aa(1):bb(1)) = ...
            cmprob(pp,1);
        
        spiim(aa(2):bb(2), aa(1):bb(1)) = ...
            spimat(pp,1);
        c1im(aa(2):bb(2), aa(1):bb(1)) = ...
            c1(pp,1);
        c2im(aa(2):bb(2), aa(1):bb(1)) = ...
            c2(pp,1);
        fuim(aa(2):bb(2), aa(1):bb(1)) = ...
            fusemat(pp,1);
        spim(aa(2):bb(2), aa(1):bb(1)) = ...
            compness(pp,1);
    end
    
    figure(1), surf( acim,'FaceLighting','phong','FaceColor','interp',...
    'EdgeColor','none');
    figure(2), surf( bcim,'FaceLighting','phong','FaceColor','interp',...
    'EdgeColor','none');
    figure(3), surf( lcim,'FaceLighting','phong','FaceColor','interp',...
    'EdgeColor','none');
    figure(4), surf( fuim,'FaceLighting','phong','FaceColor','interp',...
    'EdgeColor','none');
    figure(5), surf( c1im,'FaceLighting','phong','FaceColor','interp',...
    'EdgeColor','none');
    figure(6), surf( c2im,'FaceLighting','phong','FaceColor','interp',...
    'EdgeColor','none');
    figure(10), surf( gwim,'FaceLighting','phong','FaceColor','interp',...
    'EdgeColor','none');
    figure(7), surf( spim,'FaceLighting','phong','FaceColor','interp',...
    'EdgeColor','none');
    figure(9),surf( spiim,'FaceLighting','phong','FaceColor','interp',...
    'EdgeColor','none');
end
    salm=mat2gray(spim);
if verbose >=1    
    
    cim=ori_im;
    cim(:,:,1)=salm;
    cim(:,:,2)=salm;
    cim(:,:,3)=salm;
    
    dim=ori_im;
    dim(:,:,1)=salm;
    
    figure(8), imshow([ori_im, cim, dim]);
    pause(0.001);
end
    salmap(ii,1)={salm};
end

end

