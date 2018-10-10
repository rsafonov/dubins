function [] = dubins_const()

    global LSL LSR RSL RSR RLR LRL TSL TSR LST RST TST LSEG SSEG RSEG TSEG INF TYPENAME DIRNAME DIRDATA;
    
    LSL = 0;
    LSR = 1;
    RSL = 2;
    RSR = 3;
    RLR = 4;
    LRL = 5;
    
    TSL = 6;
    TSR = 7;
    LST = 8;
    RST = 9;
    TST = 10;
    
    LSEG = 0;
    SSEG = 1;
    RSEG = 2;
    TSEG = 3;
    
    INF=100000.0;
    
    TYPENAME = ['LSL'; 'LSR'; 'RSL'; 'RSR'; 'RLR'; 'LRL'; 'TSL'; 'TSR'; 'LST'; 'RST'; 'TST'];
    DIRNAME = ['LSEG'; 'SSEG'; 'RSEG'; 'TSEG'];
    
    DIRDATA = [LSEG, SSEG, LSEG; LSEG, SSEG, RSEG; RSEG, SSEG, LSEG; RSEG, SSEG, RSEG; RSEG, LSEG, RSEG; LSEG, RSEG, LSEG; ...
        TSEG, SSEG, LSEG; TSEG, SSEG, RSEG; LSEG, SSEG, TSEG; RSEG, SSEG, TSEG; TSEG, SSEG, TSEG;];
end