import numpy as np
from wmpl.Utils.TrajConversions import latLonAlt2ECEF
from wmpl.Utils.Math import vectNorm


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
    if az < 0:
        az += 360

    return az, el, distance


if __name__ == "__main__":
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
