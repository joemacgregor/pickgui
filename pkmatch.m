function                    pkmatch
% PKMATCH Match layer numbers of picks files without having to reload the files.
% 
%   PKMATCH is an interactive function that corrects layer numbers in a
%   transect's picks files subsequent to a chosen picks file.
% 
% Joe MacGregor (UTIG)
% Last updated: 02/26/13

[file_pk, path_pk]          = uigetfile('*.mat', 'Choose last picks file of transect that contains new picks:');
if ~file_pk
    disp('No picks file selected. PKMATCH cancelled.')
    return
end

file_all                    = dir([path_pk '*.mat']);
file_all                    = {file_all.name};
for ii = 1:length(file_all)
    if strcmp(file_pk, file_all{ii})
        break
    end
end
ind_file_curr               = ii;

if (ind_file_curr == length(file_all))
    disp(['Selected picks file (' file_pk ') is the last for this transect. No later picks to adjust. PKMATCH cancelled.'])
    return
end

try
    pk_ref                  = load([path_pk file_pk]);
    pk_ref                  = pk_ref.pk;
catch %#ok<CTCH>
    disp('Something went wrong loading selected picks file. PKMATCH cancelled.')
    return
end

for ii = (ind_file_curr + 1):length(file_all)
    try
        disp(['Updating ' file_all{ii} '...'])
        pk                  = load([path_pk file_all{ii}]);
        pk                  = pk.pk;
        pk.ind_match        = NaN(pk.num_layer, 1);
        tmp1                = NaN(pk_ref.num_layer, pk.ind_overlap(1));
        for jj = 1:pk_ref.num_layer
            tmp1(jj, :)     = pk_ref.layer(jj).twtt_smooth(pk_ref.ind_overlap(2):end);
        end
        for jj = 1:pk.num_layer
            tmp2            = NaN(pk_ref.num_layer, pk.ind_overlap(1));
            tmp2(:, ~isnan(pk.layer(jj).twtt_smooth(1:pk.ind_overlap(1)))) ...
                            = tmp1(:, ~isnan(pk.layer(jj).twtt_smooth(1:pk.ind_overlap(1)))) - repmat(pk.layer(jj).twtt_smooth(~isnan(pk.layer(jj).twtt_smooth(1:pk.ind_overlap(1)))), pk_ref.num_layer, 1);
            tmp3            = NaN(1, pk_ref.num_layer);
            for kk = 1:pk_ref.num_layer
                tmp3(kk)    = abs(nanmean(tmp2(kk, :)));
            end
            if (length(find(tmp3 < pk.twtt_match)) == 1) % found a single matching layer within the threshold
                pk.ind_match(jj) ...
                            = pk_ref.ind_match(logical(tmp3 < pk.twtt_match));
            end
        end
        if any(isnan(pk.ind_match))
            pk.ind_match(isnan(pk.ind_match)) ...
                            = (max(pk_ref.ind_match) + 1):(max(pk_ref.ind_match) + length(find(isnan(pk.ind_match))));
        end        
        pk                  = orderfields(pk);
        save([path_pk, file_all{ii}], '-v7.3', 'pk')
        pk_ref              = pk;
    catch %#ok<CTCH>
        disp(['Something went wrong updating ' file_all{ii} '. PKMATCH cancelled.'])
        return
    end
end

disp(['PKMATCH finished updating all picks files in ' path_pk ' after ' file_pk '.'])