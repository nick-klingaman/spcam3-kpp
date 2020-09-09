begin
  hyai= \
  (/0.00225523952394724, 0.00503169186413288, 0.0101579474285245,  \
    0.0185553170740604, 0.0306691229343414, 0.0458674766123295,  \
    0.0633234828710556, 0.0807014182209969, 0.0949410423636436,  \
    0.11169321089983, 0.131401270627975, 0.154586806893349,  \
    0.181863352656364, 0.17459799349308, 0.166050657629967,  \
    0.155995160341263, 0.14416541159153, 0.130248308181763,  \
    0.113875567913055, 0.0946138575673103, 0.0753444507718086,  \
    0.0576589405536652, 0.0427346378564835, 0.0316426791250706,  \
    0.0252212174236774, 0.0191967375576496, 0.0136180268600583,  \
    0.00853108894079924, 0.00397881818935275, 0, 0 /);
  hybi= \
  (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0393548272550106,  \
    0.0856537595391273, 0.140122056007385, 0.204201176762581,  \
    0.279586911201477, 0.368274360895157, 0.47261056303978,  \
    0.576988518238068, 0.672786951065063, 0.753628432750702,  \
    0.813710987567902, 0.848494648933411, 0.881127893924713,  \
    0.911346435546875, 0.938901245594025, 0.963559806346893,  \
    0.985112190246582, 1 /);

ps = new((/1,1/),float)
p0 = 1000.
ps = p0

pp = pres_hybrid_ccm(ps,p0,hyai,hybi)
p = pp(:,0,0)

print(p)

pcrm = (/3.2499485,  6.35785437,  12.357132,  21.9743328,  35.3368607,  51.732933,  69.5498428,  \
86.3000565,  101.136139,  116.4142, 134.360489,  155.740967,  181.203476,  211.606659,  248.460571,  \
292.412933,  344.863068,  407.171509,  480.695465,  563.87561, 649.06665,  727.82666,  787.623535,  \
822.654297,  842.52832,  855.141357,  865.338623,  875.632874,  884.720947,  \
891.258423,  896.484802, 901.733032,  \
907.003174,  912.295227,  917.609314,  922.945496,  928.303772,  933.684265,  939.086975,  944.512024, \
949.959351,  955.429077,  960.921265,  966.435974,  971.973267,  977.533142,  983.115662,  988.721008, \
994.34906,  1000./)

nz = dimsizes(pcrm)

hyacrm = pcrm
hybcrm = pcrm

do k=0,nz-1
 do m =1,dimsizes(p)-1
   if(p(m).ge.pcrm(k)) then
    hyacrm(k) = hyai(m) + (hyai(m-1)-hyai(m))/(p(m-1)-p(m))*(pcrm(k)-p(m))
    hybcrm(k) = hybi(m) + (hybi(m-1)-hybi(m))/(p(m-1)-p(m))*(pcrm(k)-p(m))
    break
   end if
 end do
end do

pp1 = pres_hybrid_ccm(ps,p0,hyacrm,hybcrm)
p1 = pp1(:,0,0)
print(hyacrm+"         "+hybcrm+"         "+pcrm+"           "+p1)
asciiwrite("hya.48",hyacrm)
asciiwrite("hyb.48",hybcrm)

end

