function write_nc_struct(nc_file, dimensions, data, attributes, is_hdf);

if nargin < 5
  is_hdf = 0;
end

if ~is_hdf
  disp(['Writing ' nc_file]);
  ncid = netcdf.create(nc_file, netcdf.getConstant('CLOBBER'));
else
  disp(['Writing ' nc_file ' (NetCDF4/HDF5)']);
  %ncid = netcdf.create(nc_file, bitor(netcdf.getConstant('CLOBBER'),netcdf.getConstant('NETCDF4')));
  unix(['rm -f ' nc_file]);
  ncid = netcdf.create(nc_file, netcdf.getConstant('NETCDF4'));
end

% Define dimensions
dimension_names = fieldnames(dimensions);
for ii = 1:length(dimension_names)
  disp(['   Adding dimension ' dimension_names{ii}]);
  netcdf.defDim(ncid,dimension_names{ii},getfield(dimensions, dimension_names{ii}));
end

% Set global attributes
if isfield(attributes,'global')
  attribute_names = fieldnames(attributes.global);
  GLOBAL = netcdf.getConstant('GLOBAL');
  for ii = 1:length(attribute_names)
    disp(['   Adding global attribute ' attribute_names{ii}]);
    value = getfield(attributes.global, attribute_names{ii});
    netcdf.putAtt(ncid, GLOBAL, attribute_names{ii}, convert_type(value));
  end
else
  warning('No global attributes found');
end

% Set variable information
variable_names = fieldnames(data);
for ii = 1:length(variable_names)
  disp(['   Adding variable ' variable_names{ii}]);
  value = getfield(data, variable_names{ii});
  varattributes = getfield(attributes, variable_names{ii});
  vardimid = [];
  if isfield(varattributes,'dimensions')
    vardimensions = getfield(varattributes, 'dimensions');
    for jj = 1:length(vardimensions)
      vardimid(jj) = netcdf.inqDimID(ncid,vardimensions{jj});
    end
  end
  [nc_class, straightjacket] = class2ncclass(value);
  varid = netcdf.defVar(ncid,variable_names{ii},nc_class,vardimid);

  % Set attributes
  attribute_names = fieldnames(varattributes);
  for jj = 1:length(attribute_names)
    if strcmp(attribute_names{jj}, 'dimensions')
      continue;
    elseif strcmp(attribute_names{jj}, 'missing_value')
      %netcdf.defVarFill(ncid,varid,false,varattributes.missing_value);
      %disp(['      setting _FillValue = ' num2str(varattributes.missing_value)]);  
    elseif strcmp(attribute_names{jj}, 'FillValue_')
      %netcdf.defVarFill(ncid,varid,false,varattributes.FillValue_);
      %disp(['      setting _FillValue = ' num2str(varattributes.FillValue_)]);
      value = getfield(varattributes, attribute_names{jj});
      disp(['      adding _FillValue = ' num2str(value)]);
      if is_hdf
	netcdf.defVarFill(ncid, varid, false, convert_type(value));
      else
	netcdf.putAtt(ncid, varid, '_FillValue', convert_type(value));
      end
      continue;
    elseif strcmp(attribute_names{jj}, 'ChunkSizes_')
      value = getfield(varattributes, attribute_names{jj});
      netcdf.defVarChunking(ncid, varid,'CHUNKED',value);
      continue;
    elseif strcmp(attribute_names{jj}, 'DeflateLevel_')
      value = getfield(varattributes, attribute_names{jj});
      netcdf.defVarDeflate(ncid, varid, true, true, value);
      continue;
    end
    disp(['      adding attribute ' variable_names{ii} '.' attribute_names{jj}]);
    value = getfield(varattributes, attribute_names{jj});
    netcdf.putAtt(ncid, varid, attribute_names{jj}, convert_type(value));
  end
end

netcdf.endDef(ncid);

% Fill variables
variable_names = fieldnames(data);
for ii = 1:length(variable_names)
  disp(['   Filling variable ' variable_names{ii}]);
  varid = netcdf.inqVarID(ncid,variable_names{ii});
  value = getfield(data, variable_names{ii});
  varattributes = getfield(attributes, variable_names{ii});
  nc_class = class2ncclass(value);
  if isfield(varattributes,'missing_value')
    missing_value = varattributes.missing_value;
    value(find(isnan(value))) = missing_value;
  elseif isfield(varattributes,'FillValue_')
    missing_value = varattributes.FillValue_;
    value(find(isnan(value))) = missing_value;
  end
  netcdf.putVar(ncid,varid,value);
end

netcdf.close(ncid)

function [nc_class, straightjacket]  = class2ncclass(variable)

matclass = class(variable);
straightjacket = '';
if strcmp(matclass, 'double') | strcmp(matclass, 'single')
  nc_class = 'NC_FLOAT'; % Only ever output floats...
elseif strcmp(matclass, 'char')
  nc_class = 'NC_CHAR';
elseif strcmp(matclass, 'int8')
  nc_class = 'NC_BYTE';
elseif strcmp(matclass, 'int16')
  nc_class = 'NC_SHORT';
  straightjacket = 'double';
elseif strcmp(matclass, 'int32')
  nc_class = 'NC_INT';
else
  warning(['Unable to save matlab variables of type ' matclass ' in NetCDF']);
  nc_class = '';
end


function output = convert_type(variable)

matclass = class(variable);
if strcmp(matclass, 'double')
  output = single(variable); % Only ever output floats...
elseif strcmp(matclass, 'char') | strcmp(matclass, 'int8') ...
       | strcmp(matclass, 'int16') | strcmp(matclass, 'int32') | strcmp(matclass, 'single')
  output = variable;
else
  error(['Unable to save matlab variables of type ' matclass ' in NetCDF']);
end

