function SHS = SHSmap(maps)
    SHS = sum(conj(maps).*maps,3);
return;