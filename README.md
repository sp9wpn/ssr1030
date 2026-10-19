# ssr1030

Live analyser for SSR / Mode S interrogations on 1030 MHz, received with a HackRF.
It is based on [dump1030](https://github.com/Fuziih/dump1030). dump1030 only pattern-matches pulses at 2.5 MS/s.
ssr1030 samples at 16 MS/s and fully decodes the Mode S uplink DPSK.

## Build

```
sudo apt install libhackrf-dev     # build header; runtime libhackrf0 is usually present
make
```

## Run

```
./ssr1030                          # live dashboard (HackRF, 1030 MHz, 16 MS/s)
./ssr1030 -o log.csv               # + every interrogation to CSV
./ssr1030 -o log.txt --space       # + space-delimited instead
./ssr1030 -w capture.iq            # + raw IQ recording (32 MB/s) for replay
./ssr1030 -f capture.iq            # replay a recording, print summary
./ssr1030 -f capture.iq -d         # replay with the dashboard
./ssr1030 -L                       # scrolling one line per interrogation
./ssr1030 -r SI17                  # only radar SI17
./ssr1030 -r SI17 -S               # only radar SI17, Mode S interrogations only
./ssr1030 -B                       # learn callsigns from ADS-B (1090) windows (kept in memory)
./ssr1030 -r SI17 -R               # only roll-calls of SI17 (no all-calls)
./ssr1030 -A 48AF01                # only interrogations addressed to one aircraft
./ssr1030 -A LOT3KM -B             # same, by callsign
./ssr1030 --adsblol 52.23,21.01    # preload callsigns of aircraft within 150 NM from adsb.lol
```

## Filters

These apply to statistics, bursts, the log and `--list`. TCAS stays a separate statistic.

- `-r, --radar SI17|II3`: keeps one radar.
  - Mode S interrogations carry the code: UF11 has II/SI, and UF4/5/20/21 have SIS (DI=3) or IIS (DI=0/1/7).
  - A/C and all-call interrogations are attributed through the main-beam burst they belong to. A burst gets the majority code of the interrogations inside it.
  - Outside bursts, an all-call is attributed if a UF11 was sent less than 300 µs before it.
  - Anything left unattributed is dropped.
- `-S, --modes-only`: drops plain Mode A/C and intermode (A/C-only and A/C/S all-calls).
- `-R, --rollcall-only`: keeps only interrogations addressed to an aircraft (UF4/5/20/21/24/19). It drops UF11 `S`, all intermode all-calls and plain Mode A/C.
- `-A, --aircraft ID`: keeps only interrogations addressed to one aircraft (UF0/4/5/16/19/20/21/24).
  - `ID` is 6 hex digits or a callsign. Use `cs:` to force a callsign that looks like hex.
  - A callsign resolves once the ADS-B table knows it.

Burst detection still uses all interrogations, because the beam is physical.
Per-burst counts, sequences, dwell and rotation are computed only from interrogations that pass the filters.

## ADS-B callsigns

`-B, --adsb` retunes the HackRF to 1090 MHz for `--adsb-dwell` seconds (default 2) at most every `--adsb-every` seconds (default 30).

- **Timing:** once the rotation is known, a window starts right after a main-beam pass and ends before the next one, so no burst is lost. A burst that does touch a window is marked `*` and excluded from per-burst statistics.
- **Sources:**
  - DF17/18 identification (CRC must be zero).
  - DF11 replies (address presence only).
  - DF20/21 Comm-B BDS 2,0, accepted only for addresses already known.
- **Preload:** `--adsblol LAT,LON` fills the table at startup with every aircraft within `--adsblol-radius` NM (default 150, API max 250) of that position, from the public [adsb.lol](https://api.adsb.lol) API (`/v2/lat/…/lon/…/dist/…`, via `curl`, no key needed). It sends the coordinates to adsb.lol and works with or without `-B` and in `--file` mode. If the request fails, the run continues without it.
- **Table:** kept in memory only. Optionally, `--adsb-db PATH` loads it at start and saves it after each window (for example to a file in `/tmp`). That also lets `--file` replays show callsigns.

Dashboard keys: `q` quit, `+`/`-` VGA ±2 dB, `L`/`l` LNA ±8 dB, `a` amp on/off, `r` reset statistics.

Military Mode 1/2 and Mode B/D are not detected. On air they were almost all coincidences of unrelated pulses.

Gain: the default is RF amp off, LNA 32, VGA 24. At 15 km from the radar the amp was needed earlier; switch it on with `-a 1` or the `a` key if the dashboard says `WEAK SIGNAL`.
The dashboard shows the current HackRF settings (AMP, LNA, VGA, nominal gain, baseband filter).
The level check covers the last max(20 s, 3 rotations), so every warning spans several beam passes:
- `CLIPPING n s ago`: any clipping in that window.
- `WEAK SIGNAL`: the strongest pulse is less than 25 dB above noise.
- `level OK`: neither.
Changing the gain restarts the level check.
If the dashboard shows `CLIPPING`, lower the gain.
Close to the radar (a few km), start with `-l 16 -g 20`.

## What it detects

| Kind | Timing from P1 |
|---|---|
| Mode A / C | P3 at 8 / 21 µs |
| A/C-only all-call (A or C) | P4 0.8 µs, 2 µs after P3 |
| A/C/S all-call (A or C) | P4 1.6 µs |
| Mode S | P2 at 2 µs, P6 at 3.5 µs, then DPSK decode of UF0/4/5/11/16/19/20/21/24 |

Each interrogation also carries these flags:
- **SLS:** P2 ≥ P1 (A/C), or no sync phase reversal (Mode S: P5 present). Both mean a sidelobe.
- **P2 level:** P2 relative to P1, in dB.
- **S1:** a whisper-shout pulse is present.

Mode S fields decoded:
- UF4/5/20/21: PC, RR, DI, IIS/SIS, and the requested BDS (RR ≥ 16 → BDS1, RRS → BDS2).
- UF11: II/SI code and PR.
- The interrogated aircraft address.

### Address recovery

Uplink AP = parity XOR high24(A(x)·G(x)) (Gertz, MIT LL ATC-117).
So A = quotient of (parity XOR AP)·x²⁴ / G.
This was checked on air: UF11 and UF16 broadcasts decode to `FFFFFF`, and roll-calls decode to real ICAO addresses.

### Bit timing

Data phase reversals start at sync + 0.5 µs (N ≥ 2 chips, Annex 10 Vol IV 3.1.2.1.5.2.3), not at sync + 0.25 µs.

## TCAS / ACAS

UF0, UF16 and whisper-shout (S1) interrogations are counted as a receiver health statistic.
They are left out of all radar analysis.

## Bursts (antenna main beam)

- **Background:** the median of the strongest radar interrogation in each 10 ms bin, over the last 20 s. This approximates the sidelobe level.
- **Burst:** a run of interrogations ≥ `--member-db` (3 dB) above background, with gaps < `--burst-gap` (30 ms). It must contain at least one interrogation ≥ `--core-db` (6 dB) above background.
- **Radar identity:** taken from the UF11 II/SI code in the burst, or from SIS.
- **Per radar:** rotation period (median time between burst peaks), dwell time, and interrogations per burst by query type (avg/min/max).
- **Burst table columns** (each interrogation counted once, so they add up to `n`):
  - `AA/CC`, `A/C`, `ACSall`: plain Mode A/C, A/C-only all-calls, A/C/S all-calls.
  - `S`: UF11 all-calls.
  - `alt`, `sqk`: UF4 and UF5 roll-calls without a Comm-B request.
  - `acid`, `setAlt`, `trk`, `hdgSpd`: Comm-B requests for BDS 2,0, 4,0, 5,0 and 6,0, whatever the format (UF4/5/20/21).
  - `oBDS`: any other BDS.
  - `msg`: Comm-A without a request.
  - `other`: UF24, UF19 and the rest.
- **Interleaving:** the most frequent in-beam transitions, with Δt, plus the full token sequence of the last burst.
  - Tokens (`#N` = repeated N times):

| Token | Meaning |
|---|---|
| `AA` `CC` | Plain Mode A / Mode C (no P4) |
| `A` `C` | A/C-only all-call (short P4) |
| `As` `Cs` | A/C/S all-call (long P4) |
| `S` | UF11 Mode S all-call (acquisition) |
| `alt` / `sqk` | UF4 altitude / UF5 Mode A code roll-call, no Comm-B |
| `acid` `setAlt` `trk` `hdgSpd` | Roll-call requesting Comm-B BDS 2,0 aircraft identification (callsign) / 4,0 selected altitude / 5,0 track and turn / 6,0 heading and speed |
| `cap` `gicb` `pos` `met` `ra` `aicb` … | BDS 1,0 / 1,7 / 5,1-5,2 / 4,4 / 3,0 / air-initiated; `S-bXY` = other BDS X,Y |
| `S-msg` | Comm-A UF20/21 (`S-msg/hdgSpd` = also requests BDS 6,0) |
| `S-elm` `S-mil` `S-sls` | UF24 ELM, UF19, Mode S without sync (sidelobe) |

Roll-calls carry the number of their target, for example `alt/1` or `setAlt/2#2`.
Targets are listed under the sequence, numbered in alphabetical order of their callsign or hex: `targets: 1/452124 2/DLH123`.

The interleaving statistics use the same tokens, without the target numbers.

## Log columns

`t_s utc kind uf amp_dbfs p2_rel_db p4_us s1 sidelobe tcas burst radar addr callsign ii_si pr rr di bds sync quality hex`

- `burst` is -1 outside a main-beam burst.
- `radar` is the attributed radar code, or `-` when unknown.
- `quality` is the weakest DPSK bit confidence (0–1).
- `hex` is the raw Mode S message.
