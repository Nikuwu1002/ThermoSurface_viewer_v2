| Datensatz             | erzeugende Funktion       | Struktur                     | Größen        | 3D-Verwendung |
|-----------------------|---------------------------|------------------------------|---------------|---------------|
| Isothermen            | `multi_isotherm`          | Liste `(T, v, p)`            | T, v, p       | Linien bei konstantem T |
| Isochoren             | `multi_isochore`          | Liste `(T, v, p)`            | T, v, p       | Linien bei konstantem v |
| Isobaren              | `multi_isobare`           | Liste `(T, v, p)`            | T, v, p       | Linien bei konstantem p |
| Binodale roh          | `calc_saturation_curve`   | Liste `(T, p_sat, v_l, v_v)` | T, p, zwei v  | muss in zwei Äste getrennt werden |
| Flüssigkeits-Binodale | Umstrukturierung          | `(T, v_l, p_sat)`            | T,v,p         | Siedelinie |
| Dampf-Binodale        | Umstrukturierung          | `(T, v_v, p_sat)`            | T,v,p         | Taulinie |
| Spinodale             | `calc_spinodale`          | Liste `(T, v_l, p)`          | T,v,p         | Stabilitätsgrenze, p aktuell problematisch |
|                       |                           | Liste `(T, v_v, p)`          |
| kritischer Punkt      | `critical_point`          | `(T_c, p_c, v_c)`            | T,p,v         | einzelner Marker |
| Stoffdaten            | `stoffdaten`              | `(T_c, p_c, M, omega)`       | Konstanten    | Skalierung/Normierung |