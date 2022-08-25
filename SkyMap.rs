pub struct SearchResult {
    az: f64,
    alt: f64,
}
pub struct DateTime {
    year: f64,
    month: f64,
    day: f64,
    time: f64,
}
pub struct ObserverPosition {
    lat: f64,
    lng: f64,
}
pub struct StarCoordinates {
    ra: f64,
    dec: f64,
}
pub struct SkyMap {
    star_coordinates: StarCoordinates,
    observer_position: ObserverPosition,
    date_time: DateTime,
}

impl SkyMap {
    pub fn new() -> SkyMap {
        let initial_observer_position = ObserverPosition {
            lat: 0.00,
            lng: 0.00,
        };

        let initial_star_coordinates = StarCoordinates {
            ra: 0.00,
            dec: 0.00,
        };
        let initial_date_time = DateTime {
            year: 0.00,
            month: 0.00,
            day: 0.00,
            time: 0.00,
        };
        return SkyMap {
            observer_position: initial_observer_position,
            date_time: initial_date_time,
            star_coordinates: initial_star_coordinates,
        };
    }
    pub fn rad2deg(value: f64) -> f64 {
        return value * 180.00 / 3.14159265358979;
    }
    pub fn deg2rad(value: f64) -> f64 {
        return value * 3.14159265358979 / 180.00;
    }
    pub fn set_date_time(&mut self, value: DateTime) {
        self.date_time = value;
    }
    pub fn set_observer_position(&mut self, value: ObserverPosition) {
        self.observer_position = value;
    }
    pub fn set_star_coordinates(&mut self, value: StarCoordinates) {
        self.star_coordinates = value;
    }
    pub fn calculate_local_sidereal_time(j2000: f64, time: f64, longitude: f64) -> f64 {
        let mut lst = 100.46 + 0.985647 * j2000 + longitude + 15.00 * time;
        if lst < 0.00 {
            lst += 360.00;
        } else if lst > 360.00 {
            lst -= 360.00;
        }
        return lst;
    }
    pub fn calculate_j2000(datetime: &DateTime) -> f64 {
        let jd = 367.00 * datetime.year
            - f64::floor(
                7.00 * (datetime.year + f64::floor((datetime.month + 9.00) / 12.00)) / 4.00,
            )
            - f64::floor(
                3.00 * (f64::floor((datetime.year + (datetime.month - 9.00) / 7.00) / 100.00)
                    + 1.00)
                    / 4.00,
            )
            + f64::floor(275.00 * datetime.month / 9.00)
            + datetime.day
            + 1721028.5
            + datetime.time / 24.00;
        return jd - 2451545.00;
    }
    pub fn deg2h(value: f64) -> f64 {
        return value / 15.00;
    }
    pub fn h2deg(value: f64) -> f64 {
        return value * 15.00;
    }
    pub fn asind(value: f64) -> f64 {
        return f64::asin(SkyMap::rad2deg(value));
    }
    pub fn acosd(value: f64) -> f64 {
        return f64::acos(SkyMap::rad2deg(value));
    }
    pub fn calculate_hour_angle(lst: f64, ra: f64) -> f64 {
        let mut ha = lst - ra;
        if ha < 0.00 {
            ha += 360.00;
        } else if ha > 360.00 {
            ha -= 360.00;
        }
        return ha;
    }
    pub fn calculate_az_alt(ha: f64, dec: f64, lat: f64) -> SearchResult {
        /*  math behind calculations -- conversion from HA and DEC to alt and  az
        sin(alt) = sin(DEC) * sin(LAT) + cos(DEC) * cos(LAT) * cos(HA)
        alt = asin(alt)
                   sin(DEC) - sin(alt) * sin(LAT)
        cos(a) = ---------------------------------
                        cos(alt) * cos(LAT)
        a = acos(a)
        If sin(HA) is negative,then az = a, otherwise az = 360 - a */

        let sin_dec = f64::sin(SkyMap::deg2rad(dec));
        let sin_ha = f64::sin(SkyMap::deg2rad(ha));
        let sin_lat = f64::sin(SkyMap::deg2rad(lat));
        let cos_dec = f64::cos(SkyMap::deg2rad(dec));
        let cos_ha = f64::cos(SkyMap::deg2rad(ha));
        let cos_lat = f64::cos(SkyMap::deg2rad(lat));
        let sin_alt = (sin_dec * sin_lat) + (cos_dec * cos_lat * cos_ha);
        let mut alt = f64::asin(sin_alt);
        let cos_alt = f64::cos(alt);
        let cos_a = (sin_dec - sin_alt * sin_lat) / (cos_alt * cos_lat);
        let mut a = f64::acos(cos_a);
        a = SkyMap::rad2deg(a);
        alt = SkyMap::rad2deg(alt);

        let mut _az = 0.00;
        if sin_ha > 0.00 {
            _az = 360.00 - a;
        } else {
            _az = a;
        }
        return SearchResult { az: _az, alt: alt };
    }

    pub fn calculate(&mut self) -> SearchResult {
        let j2000 = SkyMap::calculate_j2000(&self.date_time);
        let lst = SkyMap::calculate_local_sidereal_time(
            j2000,
            self.date_time.time,
            self.observer_position.lng,
        );

        let ha = SkyMap::calculate_hour_angle(lst, self.star_coordinates.ra);
        return SkyMap::calculate_az_alt(ha, self.star_coordinates.dec, self.observer_position.lat);
    }
}
fn main() {
    let mut _skymap = SkyMap::new();

    let time_and_date_of_observation = DateTime {
        year: 2021.00,
        month: 9.00,
        day: 4.00,
        time: 20.2, // UTC
    };
    //Sirius
    let star_coordinates = StarCoordinates {
        ra: 101.52,
        dec: -16.7424,
    };
    // Los Angeles
    let my_position = ObserverPosition {
        lat: 34.05,
        lng: -118.24358,
    };
    _skymap.set_date_time(time_and_date_of_observation);
    _skymap.set_star_coordinates(star_coordinates);
    _skymap.set_observer_position(my_position);
    let search_result = _skymap.calculate();
    println!("Alt {:.} Az {:.}", search_result.alt, search_result.az);
}
