##Place to add functions commonly used by the EPW scripts


##------------------------------------------------------------------------------
##Match for EPW fields

get_field_index <- function(var.name) {

   field.names <- c('year', 'month', 'day', 'hour', 'minute',
      'data_source_and_uncertainty_flags', 'dry_bulb_temperature',
      'dew_point_temperature', 'relative_humidity',
      'atmospheric_station_pressure', 'extraterrestrial_horizontal_radiation',
      'extraterrestrial_direct_normal_radition',
      'horizontal_infrared_radiation_intensity', 'global_horizontal_radiation',
      'direct_normal_radiation', 'diffuse_horizontal_radiation',
      'global_horizontal_illuminance', 'direct_normal_illuminance',
      'diffuse_horizontal_illuminance', 'zenith_luminance', 'wind_direction',
      'wind_speed', 'total_sky_cover', 'opaque_sky_cover', 'visibility',
      'ceiling_height', 'present_weather_observation', 'present_weather_codes',
      'precipitable_water', 'aerosol_optical_depth', 'snow_depth',
      'days_since_last_snowfall', 'albedo', 'liquid_precipitation_depth',
      'liquid_precipitation_quantity')
   ix <- grep(var.name,field.names)
}

##------------------------------------------------------------------------------
##Return the short names for the EPW field names
   
get_short_name <- function(var.name) {

    short.name <- switch(var.name, 
        dry_bulb_temperature='TAS',
        dew_point_temperature='DWPT',
        relative_humidity='RHS',
        atmospheric_station_pressure='PS',
        extraterrestrial_horizontal_radiation='ETHR',
        extraterrestrial_direct_normal_radition='ETDR',
        horizontal_infrared_radiation_intensity='IR',
        global_horizontal_radiation='GHR',
        direct_normal_radiation='DNR',
        diffuse_horizontal_radiation='DHR',
        wind_direction='WDIR',
        wind_speed='WSPD',
        total_sky_cover='TSC',
        opaque_sky_cover='OSC',
        snow_depth='SND',
        liquid_precipitation_quantity='PR')
   return(short.name)
}


##------------------------------------------------------------------------------
##Return the lon/lat coordinates for the supplied EPW file

get_epw_coordinates <- function(epw.dir,epw.file) {

   epw <- read.epw.file(epw.dir,epw.file)
   epw.header <- epw$header
   epw.first <- strsplit(epw.header[1],',')[[1]]
   lon <- as.numeric(epw.first[8]) ##Fixed location
   lat <- as.numeric(epw.first[7])
   if (lon < -180 | lon > 180) {
      stop('Ill defined longitude coordinate')
   }
   if (lat < 40 | lat > 90) {
      stop('Ill defined latitude coordinate')
   }

   rv <- c(lon,lat)
   return(rv)
}
