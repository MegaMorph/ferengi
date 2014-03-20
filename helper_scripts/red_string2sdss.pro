;convert the ts_file string into run, camcol, rerun, field numbers
PRO red_string2sdss, ts_file, run, camcol, rerun, field
  run = ulong(strmid(ts_file, 8, 6))
  camcol = ulong(strmid(ts_file, 15, 1))
  rerun = ulong(strmid(ts_file, 17, 2))
  field = ulong(strmid(ts_file, 20, 4))
END
