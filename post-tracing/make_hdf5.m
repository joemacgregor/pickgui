% MAKE_HDF5 Make HDF5 database of Greenland radiostratigraphy v2.
% 
% Sticking with .mat for now...HDF5 formal format way too complicated.
% 
% Joe MacGregor
% Last updated: 18 November 2024

clear

dir_mat						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
dir_breakout				= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/breakout/';
file_save					= 'Greenland_radiostratigraphy_v2.mat';
save_ver					= 2.0;

load([dir_mat 'xyz_all.mat'], 'campaign')
load([dir_mat 'pk_cat.mat'], 'depth', 'dist', 'elev_wgs84', 'file_pk', 'gps_time', 'ind_campaign', 'ind_trace_layer', 'int', 'lat', 'lon', 'num_file_pk', 'num_layer', 'num_trace', 'thick', 'twtt_ice', 'x', 'y')
load([dir_mat 'date_all.mat'], 'age', 'age_uncert')

% radar system compilation (30 campaigns, ingoring ground)
radar_system				= {'ICORDS' 'ICORDS' 'ICORDS' 'ICORDS' 'ICORDSv2' 'ICORDSv2' 'ICORDSv2' 'ICORDSv2' 'ACORDS' 'ACORDS' 'MCRDS' 'MCRDS' '' 'MCRDS' 'MCRDS' 'MCoRDS' 'MCoRDS' 'MCoRDSv2' 'MCoRDSv2' 'MCoRDSv2' 'MCoRDSv3' 'MCoRDSv3' 'MCoRDSv5' '' 'MCoRDSv5' 'MCoRDSv5' '' 'MCoRDSv3' 'MCoRDSv3' 'MCoRDSv3'};

% load mat/merge_all elev_bed_gimp elev_smooth_gimp elev_surf_gimp int_smooth twtt_smooth twtt_surf

FDM							= load([dir_mat 'gsfc_fdm_121.mat'], 'FAC', 'time_dt', 'x', 'y');
FDM.num_time_dt				= length(FDM.time_dt);
[FDM.x, FDM.y]				= deal(double(FDM.x), double(FDM.y));

% simplify filenames and unsparsify depth (much simpler), and correct for firn
[depth_pk, elev_pk, FDM_corr, int_pk, twtt_pk] ...
							= deal(cell(1, num_file_pk));
for ii = 1:num_file_pk
	[depth_pk{ii}, elev_pk{ii}, int_pk{ii}, twtt_pk{ii}] ...
							= deal(NaN(num_layer(ii), num_trace(ii)));
	depth_tmp				= full(depth{ii});
	depth_tmp(~depth_tmp)	= NaN;
	depth_pk{ii}(:, ind_trace_layer{ii}) ...
							= depth_tmp;
	elev_tmp				= full(elev_wgs84{ii});
	elev_tmp(~elev_tmp)	= NaN;
	elev_pk{ii}(:, ind_trace_layer{ii}) ...
							= elev_tmp;
	int_tmp					= full(int{ii});
	if ~isempty(find(isnan(int_tmp), 1))
		int_tmp(isnan(int_tmp)) ... 
							= 0;
	end
	int_tmp(~int_tmp)		= NaN;
	int_pk{ii}(:, ind_trace_layer{ii}) ...
							= int_tmp;
	twtt_tmp				= full(twtt_ice{ii});
	twtt_tmp(~twtt_tmp)		= NaN;
	twtt_pk{ii}(:, ind_trace_layer{ii}) ...
							= twtt_tmp;
	ind_FDM_curr			= interp1(FDM.time_dt, 1:FDM.num_time_dt, datetime(str2double(file_pk{ii}(6:13)), 'ConvertFrom', 'yyyymmdd'), 'nearest'); % time index for flight in FDM time vector
	FDM_corr{ii}			= interp2(FDM.x, FDM.y, squeeze(FDM.FAC(:, :, ind_FDM_curr)), x{ii}, y{ii}, 'linear', 0);
end
clear depth elev_wgs84 int twtt_ice

ind_campaign_ignore			= [13 24 26 27]; % campaigns to ignore (not NASA or NSF)
ind_campaign_ref			= setdiff(1:length(campaign), ind_campaign_ignore);
campaign(ind_campaign_ignore) ...
							= [];
num_campaign				= length(campaign);

%% build structure

grs2						= struct;
grs2.num_campaign			= num_campaign;
grs2.campaign				= struct;

for ii = 1:num_campaign
	disp(campaign{ii})
    grs2.campaign(ii).name	= campaign{ii};
	kk						= 1;
	grs2.campaign(ii).segment=struct;	
	grs2.campaign(ii).num_segment ...
                            = length(find(ind_campaign == ind_campaign_ref(ii)));
	for jj = find(ind_campaign == ind_campaign_ref(ii))
        grs2.campaign(ii).segment(kk).name ...
                            = file_pk{jj}(6:(end - 7));
		[grs2.campaign(ii).segment(kk).dist, grs2.campaign(ii).segment(kk).lat, grs2.campaign(ii).segment(kk).lon, grs2.campaign(ii).segment(kk).num_layer, grs2.campaign(ii).segment(kk).num_trace, grs2.campaign(ii).segment(kk).thick, ...
         grs2.campaign(ii).segment(kk).time, grs2.campaign(ii).segment(kk).x, grs2.campaign(ii).segment(kk).y, grs2.campaign(ii).segment(kk).firn_corr] ...
                            = deal(dist{jj}, lat{jj}, lon{jj}, num_layer(jj), num_trace(jj), thick{jj}, gps_time{jj}, x{jj}, y{jj}, FDM_corr{jj});
		grs2.campaign(ii).segment(kk).stratigraphy ...
                            = struct;
		for ll = 1:num_layer(jj)
			[grs2.campaign(ii).segment(kk).stratigraphy(ll).age, grs2.campaign(ii).segment(kk).stratigraphy(ll).age_uncert, grs2.campaign(ii).segment(kk).stratigraphy(ll).depth, ...
             grs2.campaign(ii).segment(kk).stratigraphy(ll).int, grs2.campaign(ii).segment(kk).stratigraphy(ll).elev, grs2.campaign(ii).segment(kk).stratigraphy(ll).twtt] ...
                            = deal((1e-3 * age{jj}(ll)), (1e-3 .* age_uncert{jj}(ll)), depth_pk{jj}(ll, :), int_pk{jj}(ll, :), elev_pk{jj}(ll, :), twtt_pk{jj}(ll, :));
            grs2.campaign(ii).segment(kk).stratigraphy ...
                            = orderfields(grs2.campaign(ii).segment(kk).stratigraphy);
		end
		grs2.campaign(ii).segment(kk) ...
                            = orderfields(grs2.campaign(ii).segment(kk));
		kk					= kk + 1;
	end
end

% global attributes
grs2.title					= 'Greenland Ice Sheet radiostratigraphy version 2 from VHF radar-sounding data collected by The University of Kansas between 1993 and 2019';
grs2.filename				= file_save;
grs2.version				= save_ver;
grs2.citation				= 'MacGregor et al., in prep., A revised and expanded deep radiostratigraphy of the Greenland Ice Sheet from airborne radar sounding surveys between 1993–2019';
grs2.date_generated			= char(datetime('now'));
grs2.point_of_contact		= 'Joseph MacGregor, joseph.a.macgregor@nasa.gov';
grs2.projection				= 'WGS 84 / NSIDC Sea Ice Polar Stereographic North (EPSG:3413)';
%							    NAME				UNITS	DESCRIPTION
vars2log					= {'age'                'ka'    'reflection age';
                               'age_uncert'         'ka'    'reflection age uncertainty';
                               'dist'				'm'		'distance along segment';
                               'depth'              'm'     'ice-equivalent depth to reflection, not firn-corrected)';
                               'int'				'dB'    'reflection relative echo intensity';
                               'elev'				'm'     'reflection elevation referenced to WGS84 ellipsoid';
							   'firn_corr'			'm'		'firn correction from GSFC FDM'
                               'lat'				'deg'   'latitude';
                               'lon'				'deg'   'longitude';
                               'time'				's'     'GPS measurement time since 1 January 0000';
                               'thick'				'm'     'ice thickness';
                               'twtt'				's'     'englacial traveltime to reflection';
                               'x'                  'm'		'x, EPSG:3413';
                               'y'                  'm'		'y, EPSG:3413'};
grs2.variables				= struct;
for ii = 1:size(vars2log, 1)
	grs2.variables(ii).name	= vars2log(ii, 1);
	grs2.variables(ii).units= vars2log(ii, 2);
	grs2.variables(ii).description ...
							= vars2log(ii, 3);
end

grs2						= orderfields(grs2);
save([dir_save file_save], 'grs2')

%% breakout individual campaign files (otherwise mat73.loadmat in Python complains)

vars2copy					= {'citation' 'date_generated' 'filename' 'point_of_contact' 'projection' 'title' 'variables' 'version'};
for ii = 1:length(vars2copy)
	eval([vars2copy{ii} ' = grs2.' vars2copy{ii} ';'])
end

for ii = 1:grs2.num_campaign
	disp(ii)
	campaign				= grs2.campaign(ii).name;
	segment					= grs2.campaign(ii).segment;
	num_segment				= grs2.campaign(ii).num_segment;
	save([dir_breakout file_save(1:(end - 4)) '_' grs2.campaign(ii).name '.mat'], 'citation', 'date_generated', 'filename', 'point_of_contact', 'projection', 'title', 'variables', 'version', 'campaign', 'segment', 'num_segment')
end