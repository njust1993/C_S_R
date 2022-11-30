function    a = loaddata(fpath, ds)
      fname = [fpath, ds, '.mat'];
      load(fname);
      a = Problem.A;
end

