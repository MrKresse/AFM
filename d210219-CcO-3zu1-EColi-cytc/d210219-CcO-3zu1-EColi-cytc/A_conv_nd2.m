clear;

input_dir='./incoming/';

dir_cont=dir(strcat(input_dir,'*.nd2'));

for iF=1:length(dir_cont)
% for iF=1:1
    base_str=dir_cont(iF).name;
    tmp=imreadBFmeta(strcat(input_dir,base_str));
    xdim=tmp.width;
    ydim=tmp.height;
    Ndim=tmp.zsize;
    if Ndim>1
        use_zsize=1;
    else
        Ndim=tmp.nframes;
        use_zsize=0;
    end
       
    for iN=1:Ndim
        display(strcat('detecting: ',num2str(iN),', file: ',base_str));
        
        if mod(iN,2000)==1
            clear tmp;
            if iN+2000-1<Ndim
                if use_zsize==1
                    tmp=imreadBF(strcat(input_dir,base_str),(iN:iN+2000-1),1,1);
                else
                    tmp=imreadBF(strcat(input_dir,base_str),1,(iN:iN+2000-1),1);
                end
            else
                if use_zsize==1
                    tmp=imreadBF(strcat(input_dir,base_str),(iN:Ndim),1,1);
                else
                    tmp=imreadBF(strcat(input_dir,base_str),1,(iN:Ndim),1);
                end
            end
%             tmp_back=mean(tmp,3);
            file_ext=iN;
        end
        
        if mod(iN,2000)==0
            tmp_slice=tmp(:,:,2000);
        else
            tmp_slice=tmp(:,:,mod(iN,2000));
        end
        tmp_slice=reshape(tmp_slice,xdim,ydim);
%         tmp_slice=tmp_slice-tmp_back;
               
        if iN==1
            imwrite(uint16(reshape(double(tmp_slice),ydim,xdim)),strcat('./incoming/',base_str,'.',num2str(file_ext),'.tif'),'WriteMode','overwrite');
        else
            imwrite(uint16(reshape(double(tmp_slice),ydim,xdim)),strcat('./incoming/',base_str,'.',num2str(file_ext),'.tif'),'WriteMode','append');
        end
    end
end
