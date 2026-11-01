# BDS 1,0: data link capability report

The aircraft's Mode S transponder fills this register to say what it can do. ssr1030 decodes it when a radar asks for BDS 1,0. bds_logger decodes it from any Comm-B reply, because the register names itself in its first byte.

ssr1030 and bds_logger print the fields in the groups below:

```
BDS1,0  OVC subnet v5 SSS  ELM up 0 down 0  ident SCS SI 1,7-chg  ACAS RA DO-185A  DTE 0000
```

- **MB:** the bit numbers within the 56-bit Comm-B data field. MB bit k is message bit k + 32.
- **Shown as:** the label printed when the bit is set, or the value.
- **EPKT:** the share of the 281 aircraft that reported 1,0 in a 9.5 h recording near EPKT, one reply per aircraft.

Sources: ICAO Doc 9871 Table A-2-16 as given by pyModeS, and the Eurocontrol "Mode-S Specific Services" tables.

## Header

| MB | Field | Meaning | Shown as | EPKT |
|---|---|---|---|---|
| 1–8 | BDS code | Always `0001 0000` (0x10). The register names itself, so a 1,0 reply is recognised without knowing what the radar asked for. | (checked, not shown) | all |
| 10–14 | Reserved | Must be 0. bds_logger rejects a reply otherwise. | — | — |

## Transponder

| MB | Field | Meaning | Shown as | EPKT |
|---|---|---|---|---|
| 15 | Overlay command capability (OVC) | Supports data parity: it XORs the requested register into the parity (the `DP` in ssr1030's parity column). Required for transponders certified from 2020 (DO-181E). | `OVC` | 99% |
| 17–23 | Mode S subnetwork version | Which edition of the register formats the transponder follows. 0 = no subnetwork, 4 = Doc 9871 Ed 1, 5 = Doc 9871 Ed 2. The values 1–3 are older Doc 9688 editions (exact dates not verified). | `subnet vN` | v5: 278, v0: 3 |
| 24 | Enhanced protocol indicator | 1 = Level 5 transponder: it can run Comm-B and ELM sessions with several ground stations at once. 0 = Level 2–4. | `L5` | 1 aircraft |
| 25 | Mode S specific services | Supports GICB register extraction beyond the basic ones. A radar can read 4,0, 5,0, 6,0 and so on only if this is set. | `SSS` | all |

Version 5 means Doc 9871 Ed 2 formats. That favours the pyModeS (Ed 2) layout of BDS 4,4 over the older Eurocontrol table. It is a hint, not proof.

## ELM (extended length messages)

| MB | Field | Meaning | Shown as | EPKT |
|---|---|---|---|---|
| 26–28 | Uplink ELM throughput | 0 = none. 1–6 = 16 segments every 1 s, 500 ms, 250 ms, 125 ms, 60 ms or 30 ms. 7 = unassigned. | `ELM up N` | 0 for 280 |
| 29–32 | Downlink ELM throughput | 0 = none. 1–3 = one 4-, 8- or 16-segment message per second. 4–6 = 16 segments every 500, 250 or 125 ms. 7–15 = unassigned. | `down N` | 0 for 280 |

ELM was built for the Mode S data link (ATN), which was never deployed. One aircraft reported `up 1 down 3`.

## Surveillance

| MB | Field | Meaning | Shown as | EPKT |
|---|---|---|---|---|
| 33 | Aircraft identification capability | The callsign comes from an automatic interface (FMS or control panel), so BDS 2,0 holds a real flight ID. | `ident` | all |
| 34 | Squitter capability subfield (SCS) | 1 = extended squitter position registers 0,5 and 0,6 were updated within the last ~10 s. | `SCS` | 43% |
| 35 | Surveillance identifier (SI) code capability | Understands SI codes (DI=3, SIS field), not only II codes. A radar with an SI code (like SI55 or SI17) can lock out only transponders that have this. | `SI` | all |
| 36 | Common usage GICB capability report | **Not a capability.** The bit flips each time the content of BDS 1,7 changes, which tells the radar to re-read 1,7. | `1,7-chg` | 44% |

- **SCS:** every one of these aircraft sends extended squitter, so 43% is low. Transponders differ in how they set this bit. A likely cause is that airborne aircraft don't update 0,6 (surface position).
- **1,7-chg:** 44% is close to the 50% a toggling bit would give.

## ACAS (TCAS)

| MB | Field | Meaning | Shown as | EPKT |
|---|---|---|---|---|
| 16 | ACAS operational | 1 = TCAS installed and operating. 0 = absent, failed or on standby. | `ACAS` / `no ACAS` | 99% |
| 37 | Hybrid surveillance | TCAS uses the intruder's ADS-B position to replace most of its own 1030 MHz interrogations. This lowers channel load (part of ssr1030's `1030 busy`). | `HS` | 26% |
| 38 | TA and RA | 1 = TCAS issues traffic and resolution advisories (TCAS II). 0 = traffic advisories only. | `RA` | 272 of 279 |
| 39–40 | RTCA version | 0 = DO-185 (TCAS 6.04A), 1 = DO-185A (7.0), 2 = DO-185B (7.1), 3 = reserved. | `DO-185`, `DO-185A`, `DO-185B` | 279 of 279 show A |

**The version field doesn't fit reality.** TCAS 7.1 has been mandatory in Europe since 2015–2017, yet every aircraft reports 7.0. Two explanations are possible, and neither is confirmed without the spec text:
- the transponder doesn't get the version from the TCAS unit and sends a fixed value;
- the two bits are numbered the other way round, in which case `01` would mean 7.1.

The HS, RA and version flags are printed only when bit 16 (ACAS operational) is set.

## Data link (DTE)

| MB | Field | Meaning | Shown as | EPKT |
|---|---|---|---|---|
| 41–56 | DTE status | One bit per data link subaddress 0–15: which Mode S data link applications are connected. | `DTE xxxx` (hex) | 0000 for all |
| 9 | Continuation flag | 1 = BDS 1,1 holds more capability data. | `+cont` | 0 |

The data link was never deployed, so the DTE field is all zeros everywhere.
