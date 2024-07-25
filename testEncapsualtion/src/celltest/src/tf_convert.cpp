#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double UAV::M2Lat(double d)
{
  // 1米对应的纬度变化约为1/111111.0度
  return home_lat_ + d / 111111.0;
}

double UAV::M2Lon(double d)
{
  // 根据当前纬度计算1米对应的经度变化
  double latitudeInRadians = home_lat_ * M_PI / 180.0;
  double metersPerDegreeLongitude = 111320.0 * cos(latitudeInRadians);
  return home_lon_ + d / metersPerDegreeLongitude;
}

double UAV::Lat2M(double latitude)
{
    // 1米对应的纬度变化约为1/111111.0度
    return (latitude - home_lat_) * 111111.0;
}

double UAV::Lon2M(double longitude)
{
    // 根据当前纬度计算1米对应的经度变化
    double latitudeInRadians = home_lat_ * M_PI / 180.0;
    double metersPerDegreeLongitude = 111320.0 * cos(latitudeInRadians);
    return (longitude - home_lon_) * metersPerDegreeLongitude;
}