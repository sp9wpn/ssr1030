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
./ssr1030 -C                       # only plain Mode A/C (older SSR radars)
./ssr1030 -A 48AF01                # only interrogations addressed to one aircraft
./ssr1030 -A 48AF01,LOT3KM,DLH123  # ... or to any of several
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
- `-C, --ac-only`: keeps only plain Mode A/C (tokens `AA`/`CC`, no P4), as sent by older SSR radars. It drops Mode S and all intermode all-calls (`A`/`C`, `As`/`Cs`). It cannot be combined with `-S`, `-R` or `-A`.
- `-R, --rollcall-only`: keeps only interrogations addressed to an aircraft (UF4/5/20/21/24/19). It drops UF11 `S`, all intermode all-calls and plain Mode A/C.
- `-A, --aircraft LIST`: keeps only interrogations addressed to the listed aircraft (UF0/4/5/16/19/20/21/24).
  - `LIST` holds up to 16 entries separated by commas, spaces or semicolons, for example `-A 48AF01,LOT3KM,DLH123`.
  - Each entry is 6 hex digits or a callsign. Use `cs:` to force a callsign that looks like hex.
  - A callsign resolves once the callsign table knows it. Until then the summary shows it as `LOT3KM(?)`.

Burst detection still uses all interrogations, because the beam is physical.
Per-burst counts, sequences, dwell and rotation are computed only from interrogations that pass the filters.

## ADS-B callsigns

`-B, --adsb` retunes the HackRF to 1090 MHz for `--adsb-dwell` seconds (default 2.9) at most every `--adsb-every` seconds (default 30).

- **Timing:** once the rotation is known, a window starts right after a main-beam pass and ends before the next one, so no burst is lost. A burst that does touch a window is marked `*` and excluded from per-burst statistics.
- **Sources:**
  - DF17/18 identification (CRC must be zero).
  - DF11 replies (address presence only).
  - DF20/21 Comm-B BDS 2,0, accepted only for addresses already known.
- **Preload:** `--adsblol LAT,LON` fills the table at startup with every aircraft within `--adsblol-radius` NM (default 150, API max 250) of that position, from the public [adsb.lol](https://api.adsb.lol) API (`/v2/lat/…/lon/…/dist/…`, via `curl`, no key needed). It sends the coordinates to adsb.lol and works with or without `-B` and in `--file` mode. If the request fails, the run continues without it.
- **Table:** kept in memory only. Optionally, `--adsb-db PATH` loads it at start and saves it after each window (for example to a file in `/tmp`). That also lets `--file` replays show callsigns.

Dashboard keys:

| Key | Action |
|---|---|
| `q` | Quit |
| `+` / `-` | VGA ±2 dB |
| `L` / `l` | LNA ±8 dB |
| `p` | RF amp on/off |
| `g` | Automatic gain on/off |
| `r` | Radar filter: type an SI code (`SI1`–`SI63`) or II code (`II0`–`II15`), then `Enter` |
| `a` | Aircraft filter: type hex addresses or callsigns (list separated by commas or spaces), then `Enter` |
| `S` | Mode S only on/off (same as `-S`) |
| `R` | Roll-calls only on/off (same as `-R`) |
| `b` | Burst graph screen on/off (see below) |
| `x` | Reset statistics |

- **Manual gain:** `+ - L l p` switch the gain to manual.
- **Filter input:** `r` and `a` open an input line pre-filled with the current filter.
  - An empty line removes the filter.
  - `Esc` cancels, `Backspace` edits, `Ctrl-U` clears the line.
- **Statistics:** any filter change (`r`, `a`, `S`, `R`) resets them, so counts from different filters never mix. This also works when replaying with `-f … -d`.

Military Mode 1/2 and Mode B/D are not detected. On air they were almost all coincidences of unrelated pulses.

## Gain

**Automatic gain is on by default** (`--no-auto-gain` turns it off, `g` toggles it).
- **Target:** the beam peak of the main radar (`-r`, or the radar with the most bursts) is kept near −6 dBFS.
  - Peak ≥ −2 dBFS, or clipping: the gain goes down at once.
  - Peak < −14 dBFS on two passes in a row: the gain goes up.
  - No beam pass for 10 s, or 2 rotations: the gain steps up by 4–12 dB while nothing strong is heard.
- **Timing:** changes are made only right after a beam pass, never during one. The background level is shifted by the same amount, so burst detection is not disturbed.
- **Stages:** the nominal gain (amp 11 dB + LNA + VGA) is split so that VGA stays ≥ 16 dB while LNA can move. The RF amp is used only above 80 dB, because near a radar it overloads first.
- **Manual override:** pressing `+ - L l p` switches to manual gain.
- **Starting point:** `-a/-l/-g` set the starting gain. The default is amp off, LNA 32, VGA 24.

ADS-B windows use the same gain as 1030.

The dashboard shows the HackRF settings, the mode (AUTO/MANUAL) and the last gain change with its reason.
The level check covers the last max(20 s, 3 rotations) and reports clipping, a weak signal (strongest pulse < 25 dB above noise), or OK.

Because the gain changes, `amp_dbfs` alone is not comparable over time.
The log has `gain_db` (nominal gain) and `rx_db` (= `amp_dbfs − gain_db`, a gain-independent relative level).
The burst table shows the gain each burst was received with.
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
- **SLS:** P2 ≥ P1 (A/C), or P5 ≥ P6 (Mode S). Both mean a sidelobe.
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

## Burst graph (`b`)

`b` switches the whole screen to the latest beam pass of the main radar, one row per query in time order:

```
   -262.58 ms  AA                     -18.5 ░░░░░░░           P2
    -36.42 ms  setAlt       4A2182     -9.1 ████████████P5█████████████████
    -23.35 ms  A                        1.2 ██P2█████████████████████████████████████████████████████
    +47.23 ms  AA                     -20.1 ▒▒▒
```

- **Each row:** time relative to the middle of the beam pass, query token (12 characters), target (8 characters: callsign when known, else hex address), `sls` for sidelobe queries, amplitude in dBFS, and a bar whose length is the amplitude.
- **Bar character:**
  - `█` in-beam query
  - `░` sidelobe (SLS: P2 ≥ P1, or Mode S without sync)
  - `▒` other weak query that is not SLS
- **SLS pulses:** `P2` (A/C) and `P5` (Mode S) are drawn at the column matching their level. They appear inside the bar when weaker than P1, and after the bar end when stronger, which is what suppresses the transponder.
  - P5 is measured at the sync phase reversal. P6 flips phase there but P5 does not, so with z1/z2 the samples just before and after it, P5 = |z1+z2|/2 and P6 = |z1−z2|/2. This also works when a strong P5 masks the sync (`S-sls`).
  - In Mode S, P2 is part of the preamble, so it is not marked.
- **Sidelobe Mode S is decoded.** P5 ≥ P6 means a transponder would not reply (the query is not for us), but P5 covers only the sync reversal and the data chips after it are intact.
  - The bit timing is recovered within ±half a chip of the nominal sync, and the query is shown with its real type and target plus `sls`.
  - These queries carry their radar's SI/II code, so sidelobes are attributed to the radar too.
- **Time range:**
  - Sidelobes are shown from 60 % of a rotation before the pass to 60 % of a rotation after, so neighbouring passes overlap. Before the rotation is known, ±2.4 s is used.
  - Weak queries are shown only within one dwell time (30–200 ms).
  - Rows keep being added as late sidelobes arrive.
- **Scrolling:** ↑/↓ (or `j`/`k`), PgUp/PgDn, Home/End. A new burst opens just before its first beam query. `b` returns to the dashboard.
- **Filters:** the type and aircraft filters apply. Sidelobe and weak queries are hidden by `-r` only when they were attributed to a different radar.

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
| `S-elm` `S-mil` | UF24 ELM, UF19 |
| `S-sls` | Mode S sidelobe (P5 ≥ P6) whose bits could not be recovered. Normally sidelobe Mode S is fully decoded and shown with its own token and the SLS flag |
| `S-grb` | Mode S with no sync reversal and no P5 (garbled) |

Roll-calls carry the number of their target, for example `alt/1` or `setAlt/2#2`.
Targets are listed under the sequence, numbered in alphabetical order of their callsign or hex: `targets: 1/452124 2/DLH123`.

The interleaving statistics use the same tokens, without the target numbers.

## Log columns

`t_s utc kind uf amp_dbfs gain_db rx_db p2_rel_db p5_rel_db p4_us s1 sidelobe tcas burst radar addr callsign ii_si pr rr di bds sync quality hex`

- `burst` is -1 outside a main-beam burst.
- `radar` is the attributed radar code, or `-` when unknown.
- `quality` is the weakest DPSK bit confidence (0–1).
- `hex` is the raw Mode S message.

## License

BSD 3-Clause, see [LICENSE](LICENSE).
The detection approach follows [dump1030](https://github.com/Fuziih/dump1030) (BSD 3-Clause, Copyright (c) 2022, Fuziih).
