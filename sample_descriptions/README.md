# Sample descriptions — original sample submission sheets

The spreadsheets the data owners submitted with their samples to the JGI Metabolomics (`jgi/`)
and EGSB (`egsb/`) LC-MS platforms. They describe the samples and LC-MS runs in the
Semantic Metabolite Atlas datasets and are published **as submitted** — they are not cleaned,
re-formatted or harmonized. EGSB sheets were native Google Sheets and are Google's `.xlsx` export.

Layout: `jgi/<proposal>/<sha>.xlsx.gz` and `egsb/<sha>.xlsx.gz`, one gzip-compressed workbook
per sheet, where `<sha>` is the first 16 hex digits of the SHA-256 of the workbook. `gunzip -N
<file>` restores the sheet under its original filename; check it with `sha256sum`. The
SpectraCodec message embedded in each deposited LC-MS run names its folder here.

The only edits are the removal of personal contact details — e-mail addresses are replaced with
`[redacted-email]` and phone numbers with `[redacted-phone]` (including `mailto:` links and file
properties) — and, in the EGSB sheets, the replacement of EGSB project identifiers with
`[redacted-project]`. Everything else in a workbook is unchanged; a sheet with nothing to redact
is the original file byte-for-byte. Sheets are added as each experiment's data are released.
