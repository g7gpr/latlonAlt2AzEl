import numpy as np
import matplotlib.pyplot as plt
from wmpl.Utils.TrajConversions import latLonAlt2ECEF, altAz2RADec,date2JD
from wmpl.Utils.Math import vectNorm
#from RMS.Astrometry.ApplyAstrometry import xyToRaDecPP, raDecToXYPP



def latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele):
    # Convert coordinates to ECEF
    tar_ecef = np.array(latLonAlt2ECEF(np.radians(tar_lat), np.radians(tar_lon), tar_ele))
    sta_ecef = np.array(latLonAlt2ECEF(np.radians(sta_lat), np.radians(sta_lon), sta_ele))

    # Create the ecef frame up vector
    ecef_up_vector = vectNorm(np.array([0, 0, 1]))

    # Calculate the unit vector
    target_vector = (tar_ecef - sta_ecef)
    unit_vector = vectNorm(target_vector)

    # Calculate the distance
    distance = np.sqrt((np.sum(target_vector ** 2)))

    # Calculate the local frame up vector - centre of the earth up to the station
    local_frame_up_vector = vectNorm(sta_ecef)

    # local east is perpendicular to local up and ecef z axis
    local_frame_east_vector = vectNorm(np.cross(ecef_up_vector, local_frame_up_vector))

    # local north is perpendicular to local up and east
    local_frame_north_vector = vectNorm(np.cross(local_frame_up_vector, local_frame_east_vector))

    # calculate elevation - sine of the shadow cast by the unit vector onto the local frame up vector
    el = np.degrees(np.arcsin(np.dot(unit_vector, local_frame_up_vector)))

    # calculate opposite - shadow cast by the unit vector onto the local east vector
    opposite = np.dot(unit_vector, local_frame_east_vector)

    # calculate adjacent - shadow cast by the unit vector onto the local north vector
    adjacent = np.dot(unit_vector, local_frame_north_vector)

    # calculate azimuth - opposite / adjacent
    az = np.degrees(np.arctan2(opposite, adjacent))

    # only return positive azimuth values
    #if az < 0:
    #    az += 360

    return az, el, distance


def test_cases():
    # test cases
    print("")
    print("Test case S1: Southern hemisphere meteor north and east of station")
    sta_lat = -32.3  # deg
    sta_lon = 115.8  # deg
    sta_ele = 37  # meters
    tar_lat = -31.908  # deg
    tar_lon = 116.156  # deg
    tar_ele = 81.1 * 1000
    azimuth, elevation, distance = latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele)
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("Station Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(sta_lat, sta_lon, sta_ele))
    print("Target  Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(tar_lat, tar_lon, tar_ele))
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("")
    print("Test case S2: Southern hemisphere meteor south and east of station")
    sta_lat = -32.3  # deg
    sta_lon = 115.8  # deg
    sta_ele = 37.0  # meters
    tar_lat = -32.4  # deg
    tar_lon = 116.156  # deg
    tar_ele = 81.1 * 1000
    azimuth, elevation, distance = latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele)
    print("")
    print("Station Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(sta_lat, sta_lon, sta_ele))
    print("Target  Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(tar_lat, tar_lon, tar_ele))
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("")
    print("Test case S3: Southern hemisphere meteor south and west of station")
    sta_lat = -32.3  # deg
    sta_lon = 115.8  # deg
    sta_ele = 37  # meters
    tar_lat = -32.4  # deg
    tar_lon = 115.356  # deg
    tar_ele = 81.1 * 1000
    azimuth, elevation, distance = latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele)
    print("")
    print("Station Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(sta_lat, sta_lon, sta_ele))
    print("Target  Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(tar_lat, tar_lon, tar_ele))
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("")
    print("Test case S4: Southern hemisphere meteor north and west of station")
    sta_lat = -32.3  # deg
    sta_lon = 115.8  # deg
    sta_ele = 37  # meters
    tar_lat = -31.908  # deg
    tar_lon = 115.356  # deg
    tar_ele = 81.1 * 1000
    azimuth, elevation, distance = latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele)
    print("")
    print("Station Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(sta_lat, sta_lon, sta_ele))
    print("Target  Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(tar_lat, tar_lon, tar_ele))
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("")
    print("Test case N1: Northern hemisphere meteor south and east of station")
    sta_lat = 32.3  # deg
    sta_lon = 115.8  # deg
    sta_ele = 37  # meters
    tar_lat = 31.908  # deg
    tar_lon = 116.156  # deg
    tar_ele = 81.1 * 1000
    azimuth, elevation, distance = latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele)
    print("")
    print("Station Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(sta_lat, sta_lon, sta_ele))
    print("Target  Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(tar_lat, tar_lon, tar_ele))
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("")
    print("Test case N2: Northern hemisphere meteor north and east of station")
    sta_lat = 32.3  # deg
    sta_lon = 115.8  # deg
    sta_ele = 37  # meters
    tar_lat = 32.4  # deg
    tar_lon = 116.156  # deg
    tar_ele = 81.1 * 1000
    azimuth, elevation, distance = latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele)
    print("")
    print("Station Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(sta_lat, sta_lon, sta_ele))
    print("Target  Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(tar_lat, tar_lon, tar_ele))
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("")
    print("Test case N3: Northern hemisphere meteor north and west of station")
    sta_lat = 32.3  # deg
    sta_lon = 115.8  # deg
    sta_ele = 37  # meters
    tar_lat = 32.4  # deg
    tar_lon = 115.356  # deg
    tar_ele = 81.1 * 1000
    azimuth, elevation, distance = latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele)
    print("")
    print("Station Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(sta_lat, sta_lon, sta_ele))
    print("Target  Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(tar_lat, tar_lon, tar_ele))
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("")
    print("Test case N4: Northern hemisphere meteor south and west of station")
    sta_lat = 32.3  # deg
    sta_lon = 115.8  # deg
    sta_ele = 37  # meters
    tar_lat = 31.908  # deg
    tar_lon = 115.356  # deg
    tar_ele = 81.1 * 1000
    azimuth, elevation, distance = latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele)
    print("")
    print("Station Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(sta_lat, sta_lon, sta_ele))
    print("Target  Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(tar_lat, tar_lon, tar_ele))
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("")
    print("Test case Flying Doctor: from Baldivis")
    sta_lat = -32.354356  # deg
    sta_lon = 115.805961  # deg
    sta_ele = 37  # meters
    tar_lat = -32.3912  # deg
    tar_lon = 115.8351  # deg
    tar_ele = 1.943 * 1000
    azimuth, elevation, distance = latlonaltlatlonalt2azel(sta_lat, sta_lon, sta_ele, tar_lat, tar_lon, tar_ele)
    print("")
    print("Station Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(sta_lat, sta_lon, sta_ele))
    print("Target  Lat: {:.3f} Lon: {:.3f} Ele: {:.1f}".format(tar_lat, tar_lon, tar_ele))
    print("Azimuth {:.2f}, Elevation {:.2f}, Distance {:.2f}km".format(azimuth, elevation, distance / 1000))
    print("")
    print("")


if __name__ == "__main__":



    plt.style.use('ggplot')
    plt.figure(figsize=(12,12))
    fig, ax = plt.subplots()
    az_list = []
    el_list = []
    utc_list = []
    x_list = []
    y_list = []

    adsb = open("adsb.txt")

    recalibrated_platepars = json.load("platepars_all_recalibrated.json")
    pp_ref.loadfromDict(recalibrated_platepars[1], use_flat=False)

    for line in adsb:
        #print(line)
        sta_lat = -32.354356  # deg
        sta_lon = 115.805961  # deg
        sta_ele = 37.0  # meters
        splitadsbdata = line.split(",")
        #print("Lat: " + (splitadsbdata[14]) + " Lon:" + (splitadsbdata[15]))

        camera_az = 94.29
        camera_el = 49.25
        camera_fovh = 88.85
        camera_fovv = 46.96

        az_min = camera_az - camera_fovh / 2
        az_max = camera_az + camera_fovh / 2

        el_min = camera_el - camera_fovv / 2
        el_max = camera_el + camera_fovv / 2

        if len((splitadsbdata[14])) != 0 and len((splitadsbdata[15])) != 0:
         aircraft_utc = splitadsbdata[6] + " " + splitadsbdata[7]
         aircraft_alt = float(splitadsbdata[11]) / 3.2808399
         aircraft_lat = float(splitadsbdata[14])
         aircraft_lon = float(splitadsbdata[15])

         #print("Station  Latitude:{:.4f} Longitude:{:.4f} Altitude:{:.0f}".format(sta_lat, sta_lon, sta_ele))
         az, el, dist  = latlonaltlatlonalt2azel(sta_lat,sta_lon,sta_ele,aircraft_lat,aircraft_lon,aircraft_alt)




         if az_min < az < az_max and el_min < el < el_max:


          x_val = np.cos(np.radians(el)) * np.sin(np.radians(az))
          y_val = np.sin(np.radians(el))
          print("Target   Latitude:{:.4f} Longitude:{:.4f} Altitude:{:.0f}".format(aircraft_lat, aircraft_lon,
                                                                                      aircraft_alt))
          print("UTC {} Azimuth:{:.2f} Elevation:{:.2f} Distance{:.2f}km x_val {:.2f} y_va {:.2f}".format(aircraft_utc,az, el, dist / 1000,
                                                                                                   x_val, y_val))
          az_list.append(az)
          el_list.append(el)
          x_list.append(x_val)
          y_list.append(y_val)
          utc_list.append(aircraft_utc)
          aircraft_year = aircraft_utc[0:4]
          aircraft_month = aircraft_utc[5:7]
          aircraft_day = aircraft_utc[8:10]
          aircraft_hour = int(aircraft_utc[11:13])
          aircraft_minute = aircraft_utc[14:16]
          aircraft_second = aircraft_utc[17:19]
          jd=date2JD(aircraft_year,aircraft_month,aircraft_day,aircraft_hour,aircraft_minute,aircraft_second)
          aircraft_ra, aircraft_dec = altAz2RADec(np.radians(az),np.radians(el),jd,np.radians(sta_lat),np.radians(sta_lon))
          print("Ra: {:.3f} Dec: {:.3f}".format(np.degrees(aircraft_ra) ,np.degrees(aircraft_dec)))

          image_x, image_y = raDecToXYPP(aircraft_ra,aircraft_dec,jd,platepar)









    plt.title("Aircraft plotted in Az and El")
    plt.scatter(az_list,el_list)



    print("Az min {}, Az max {}".format(az_min,az_max))
    print("El min {}, El max {}".format(el_min,el_max))

    xmin = np.cos(np.radians(0)) * np.sin(np.radians(-5-camera_fovh/2))
    xmax = np.cos(np.radians(camera_el-camera_fovv/2)) * np.sin(np.radians(camera_az+camera_fovh/2))



    ymin = np.sin(np.radians(camera_el-camera_fovv/2))
    ymax = np.sin(np.radians(camera_el+camera_fovv/2))




    ax.set_xlim(70,110)
    ax.set_ylim(30,70)


    plt.xlabel('azimuth (clockwise from north)')
    plt.ylabel('elevation (degrees above horizontal)')
    for i in range(len(utc_list)):
     plt.annotate(utc_list[i], (az_list[i], el_list[i]), fontsize=4)

    plt.savefig("flight.png")
    #test_cases()