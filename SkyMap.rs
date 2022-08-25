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
    newday: f64,
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
            newday: 0.00,
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
        /*  math behind calculations -- conversion from HA and DEC to ALT and  AZ
        sin(ALT) = sin(DEC) * sin(LAT) + cos(DEC) * cos(LAT) * cos(HA)
        ALT = asin(ALT)
                   sin(DEC) - sin(ALT) * sin(LAT)
        cos(A) = ---------------------------------
                        cos(ALT) * cos(LAT)
        A = acos(A)
        If sin(HA) is negative,then AZ = A, otherwise AZ = 360 - A */

        let sinDEC = f64::sin(SkyMap::deg2rad(dec));
        let sinHA = f64::sin(SkyMap::deg2rad(ha));
        let sinLAT = f64::sin(SkyMap::deg2rad(lat));
        let cosDEC = f64::cos(SkyMap::deg2rad(dec));
        let cosHA = f64::cos(SkyMap::deg2rad(ha));
        let cosLAT = f64::cos(SkyMap::deg2rad(lat));
        let sinALT = (sinDEC * sinLAT) + (cosDEC * cosLAT * cosHA);
        let mut ALT = f64::asin(sinALT);
        let cosALT = f64::cos((ALT));
        let cosA = (sinDEC - sinALT * sinLAT) / (cosALT * cosLAT);
        let mut A = f64::acos(cosA);
        A = SkyMap::rad2deg(A);
        ALT = SkyMap::rad2deg(ALT);

        let mut AZ = 0.00;
        if sinHA > 0.00 {
            AZ = 360.00 - A;
        } else {
            AZ = A;
        }
        return SearchResult { az: AZ, alt: ALT };
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
        year: 2022.00,
        month: 12.00,
        day: 6.00,
        time: 13.45,
    };
    let star_coordinates = StarCoordinates {
        ra: 150.00,
        dec: -122.00,
    };
    let my_position = ObserverPosition {
        lat: 50.00,
        lng: 21.00,
    };
    _skymap.set_date_time(time_and_date_of_observation);
    _skymap.set_star_coordinates(star_coordinates);
    _skymap.set_observer_position(my_position);
    let search_result = _skymap.calculate();
    println!("Alt {:.} Az {:.}", search_result.alt, search_result.az);
}
