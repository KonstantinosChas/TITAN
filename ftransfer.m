function ftransfer(diskorig,diskdest,folderorig,dataorig,datadest,filext,filename)

% checks for files in designated folder(s)
% copies files in a new folder
% deletes old files
% checks again for new files in the designated folder 
% for continuous copy/delete use with while loop in parent code

folderdest = folderorig ;


dir_orig = fullfile(diskorig,folderorig,dataorig); 
dir_main = fullfile(diskorig,folderorig);
cd(dir_main);

cd(dir_orig);  % change directory to where the original images are initially being stored
listing_dest = dir(['*.',filext]);  % list images in rawdata folder (raw images captured by EBUS)

if numel(listing_dest) > 0    % check if there are new images in the pre-defined folder (new images should mean new hydrophone data as well)   
    newfoldername = (['im_',char(datetime('now', 'Format','dd_MM_HH_mm_ss'))]);  % new folder name for each new experiment : date-time relevant
    dir_dest = fullfile(diskorig,folderdest,diskdest,newfoldername);    % path for new folder to copy images in
    copyfile(['*.',filext],dir_dest); % copy images to new folder: newfoldername
    
    % pre-allocate
    f = cell(1,length(listing_dest));
    ext = cell(1,length(listing_dest));
    old_files = cell(1,length(listing_dest));
    newfilename = cell(1,length(listing_dest));
    
    cd(dir_orig); 
    % new names of files
    for id = 1 : length(listing_dest)
        if id > 1
            if id < 10 
                newfilename{id} = ([filename,'_00',num2str(id)]); % names of images (001,002...)
            elseif id>= 10 && id<100
                newfilename{id} = ([filename,'_0',num2str(id)]);
            elseif id>= 100 && id<1000
                newfilename{id} = ([filename,'_',num2str(id)]);
            end
        else 
            newfilename{id} = [filename,'_',num2str(id)];  
        end
    end

    % rename files in new folder
    for id = 1:length(listing_dest)  % do this for all images found in the folder
        [~, f{id},ext{id}] = fileparts(listing_dest(id).name);  % get file parts
        rename = strcat(newfilename{id},ext{id}) ;  % set the new name 
        cd(dir_dest);
        movefile(listing_dest(id).name, rename);    % rename
    end

    % delete files in original folders 
    cd(dir_orig);
    listing_orig = dir(['*.',filext]); 
    for id = 1:length(listing_orig) % do this for all images found in the folder
        disp('here i am'); % indication that deleting is done 
        old_files{id} = listing_orig(id).name;   % list the old file names 
        delete(old_files{id});   % delete the old file
    end
end
end


