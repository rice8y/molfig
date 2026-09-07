//! Bounded approximation of the reference GPU's sinc-like instruction.
//!
//! The instruction discards fractional bits below Q0.23 before evaluating
//! its approximation. These 32 cubic Chebyshev segments approximate its
//! binary32 rounding cells over that entire input interval. They are not
//! indexed by a pixel, viewport, molecule, or noise seed. Evaluation remains
//! in binary64 until the instruction-result boundary.
//!
//! This is not a bit-exact SFU emulator: independently sampled inputs still
//! expose rare one-ULP instruction-result differences. In particular, do not
//! infer exact noise parity merely from the polynomial's small numeric error.

#[rustfmt::skip]
#[allow(clippy::excessive_precision)]
const COEFFICIENTS: [[f64; 4]; 32] = [
    [1.5705597703094922e+0, -3.1538187054687877e-4, -7.8838194309043986e-5, 4.2049128935021247e-9],
    [1.5692986142833758e+0, -9.4567980858462509e-4, -7.8722760421995292e-5, 1.4729804042595340e-8],
    [1.5667781312548577e+0, -1.5746150648799154e-3, -7.8495570742690630e-5, 2.4198223703735686e-8],
    [1.5630019533720743e+0, -2.2012777006050047e-3, -7.8153018192391761e-5, 3.3529882389261143e-8],
    [1.5579755458489530e+0, -2.8247515721110082e-3, -7.7698476725602598e-5, 4.2770659864351959e-8],
    [1.5517061617324690e+0, -3.4441575860736716e-3, -7.7135906905885688e-5, 5.2086619523522963e-8],
    [1.5442028607718317e+0, -4.0585865496033874e-3, -7.6457893039867685e-5, 6.1335240713750854e-8],
    [1.5354764689613105e+0, -4.6671591808478477e-3, -7.5671633320050633e-5, 7.0811332772488844e-8],
    [1.5255395662295514e+0, -5.2690022201729390e-3, -7.4774123569258964e-5, 7.8230612638696819e-8],
    [1.5144064802994419e+0, -5.8632549118081793e-3, -7.3771987962439284e-5, 8.7414108057654316e-8],
    [1.5020932372670228e+0, -6.4490661293499603e-3, -7.2665633958250597e-5, 9.6733072423417907e-8],
    [1.4886175547714779e+0, -7.0256070918809314e-3, -7.1454879905556454e-5, 1.0409604757899059e-7],
    [1.4739988078111568e+0, -7.5920465330238237e-3, -7.0143423807672207e-5, 1.1362219289536586e-7],
    [1.4582579846720551e+0, -8.1475991991265453e-3, -6.8731725071884941e-5, 1.2089626940950768e-7],
    [1.4414176497378153e+0, -8.6914710547318537e-3, -6.7222918833147589e-5, 1.3011906049192802e-7],
    [1.4235019390743899e+0, -9.2229061473477498e-3, -6.5620810899507696e-5, 1.3764102917211575e-7],
    [1.4045364562621703e+0, -9.7411608975472262e-3, -6.3929709350742320e-5, 1.4531718678937699e-7],
    [1.3845482994054910e+0, -1.0245513220900428e-2, -6.2149160984434053e-5, 1.5248611602503682e-7],
    [1.3635659518853580e+0, -1.0735275248530492e-2, -6.0279017047402949e-5, 1.5789613281818714e-7],
    [1.3416192944858565e+0, -1.1209765584007539e-2, -5.8330599142855926e-5, 1.6555955159537085e-7],
    [1.3187394854468042e+0, -1.1668346752215744e-2, -5.6303974893907077e-5, 1.7117740831917734e-7],
    [1.2949589792739404e+0, -1.2110404597375637e-2, -5.4199249128017797e-5, 1.7878543948170308e-7],
    [1.2703114081628737e+0, -1.2535346768560909e-2, -5.2023652523707959e-5, 1.8418446846194647e-7],
    [1.2448315698328856e+0, -1.2942615915495606e-2, -4.9781101591163823e-5, 1.8983531081814292e-7],
    [1.2185553428582840e+0, -1.3331683338638437e-2, -4.7474865456062569e-5, 1.9536038821599171e-7],
    [1.1915196520476641e+0, -1.3702040478572536e-2, -4.5105563754669516e-5, 1.9903802849252092e-7],
    [1.1637623571972420e+0, -1.4053229175693096e-2, -4.2684435333019904e-5, 2.0471643381897355e-7],
    [1.1353222505890568e+0, -1.4384813578899215e-2, -4.0207131839824983e-5, 2.0844797916821123e-7],
    [1.1062389541301747e+0, -1.4696385545395223e-2, -3.7681280390377624e-5, 2.1233903740021971e-7],
    [1.0765528384325151e+0, -1.4987587917383137e-2, -3.5110736551959616e-5, 2.1591405537730008e-7],
    [1.0463050082392606e+0, -1.5258070255776815e-2, -3.2503383244598868e-5, 2.1956534334992209e-7],
    [1.0155371763954106e+0, -1.5507554611073796e-2, -2.9862105066502197e-5, 2.2129517720371842e-7],
];

pub(super) fn sinc_approximation(argument: f32) -> f32 {
    debug_assert!((0.0..=1.0).contains(&argument));
    let fixed = (argument * 8_388_608.0) as u32;
    let segment = (fixed >> 18).min(31) as usize;
    let local = f64::from(fixed - ((segment as u32) << 18)) / 131_072.0 - 1.0;
    let c = COEFFICIENTS[segment];
    // Clenshaw evaluation. Keep these binary64 operations separate on every
    // target instead of allowing host-specific fused contraction.
    let twice = 2.0 * local;
    let b0 = c[1] - c[3];
    let b1 = c[2] + c[3] * twice;
    let a0 = c[0] - b1;
    let a1 = b0 + b1 * twice;
    (a0 + a1 * local) as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_fixed_point_inputs_share_the_instruction_result() {
        for fixed in [
            1,
            127,
            1024,
            65_535,
            131_072,
            1 << 21,
            1 << 22,
            (1 << 23) - 1,
        ] {
            let left = fixed as f32 / 8_388_608.0;
            let right = f32::from_bits(left.to_bits() + 1);
            assert_eq!((left * 8_388_608.0) as u32, (right * 8_388_608.0) as u32);
            assert_eq!(
                sinc_approximation(left).to_bits(),
                sinc_approximation(right).to_bits()
            );
        }
    }

    #[test]
    fn all_fixed_point_inputs_are_monotone_and_finite() {
        let mut previous = sinc_approximation(0.0);
        for fixed in 1..=8_388_608 {
            let next = sinc_approximation(fixed as f32 / 8_388_608.0);
            assert!(next <= previous);
            assert!((1.0..=std::f32::consts::FRAC_PI_2).contains(&next));
            previous = next;
        }
        assert_eq!(previous, 1.0);
    }

    #[test]
    fn approximation_stays_close_to_the_mathematical_sinc() {
        for step in 1..=65_536 {
            let argument = step as f32 / 65_536.0;
            let x = f64::from(argument);
            let exact = libm::sin(x * std::f64::consts::FRAC_PI_2) / x;
            assert!((f64::from(sinc_approximation(argument)) - exact).abs() <= 3.0e-7);
        }
    }
}
