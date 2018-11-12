%=== Test Case 1 ===%
if RunOptions.IniCondTest1 == 1;
    p0DGM = zeros(DGMMeshD.K*DGMMeshD.Np,1);
    for ii=1:DGMMeshD.K*DGMMeshD.Np;
        if norm([DGMMeshD.x(ii),DGMMeshD.y(ii)] - [0.005,0.005],2)<0.003
            p0DGM(ii) = 0.4;
        end
    end
end
%=== Test Case 2 ===%
if RunOptions.IniCondTest2 == 1;
    p0DGM = zeros(DGMMeshD.K*DGMMeshD.Np,1);
    for ii=1:DGMMeshD.K*DGMMeshD.Np;
        if norm([DGMMeshD.x(ii),DGMMeshD.y(ii)] - [0.005,0.005],2)<0.003
            p0DGM(ii) = 0.4;
        end
        if abs(DGMMeshD.x(ii) - 0.001)<0.005 && abs(DGMMeshD.y(ii) - 0.001) < 0.001;
            p0DGM(ii) = 0.4;
        end
        if abs(DGMMeshD.y(ii) - 0.001)<0.005 && abs(DGMMeshD.x(ii) - 0.001) < 0.001;
            p0DGM(ii) = 0.4;
        end
    end
end
p0DGM = reshape(p0DGM,DGMMeshD.Np,DGMMeshD.K)*RunOptions.ScalingOptical;
PLOT.InitialPressurezAxis = [0 2]; %z axis limits for plotting QPAT elastic wave propogation
PLOT.DGMPlotzAxis = [0 2]; %z axis limits for plotting QPAT elastic wave propogation
PLOT.DGMAdjointPlotzAxis = [0 2]; %z axis limits for plotting QPAT elastic wave propogation
PLOT.HoldColourAxis = 1; %Hold colour axis
PLOT.ColourAxis = [0 2]; %Colour axis