figure;

t = [0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.00]


hold on


densitysearchsgs_arr = [40.56642498130963,44.6015168647238,47.60950105311776,49.64964858286437,51.8106083088252,55.0432028507029,58.60111195707734,59.08506687854077,62.179635444985436,65.82955330321249]

plot(t,densitysearchsgs_arr, 'LineWidth', 2)





hold on

repeated_arr = [46.236684877402126,49.79459398377656,50.44772557566571,50.4922603166337,53.81594792012029,56.852886200921766,59.59722889942855,60.22408307137327,63.78317969771177,66.24133022571817]

plot(t,repeated_arr, 'LineWidth', 2)





hold on

greedy_arr = [46.236684877402126,49.79459398377656,50.44772557566571,50.4922603166337,53.81594792012029,56.852886200921766,59.59722889942855,60.22408307137327,62.31436382328442,66.24133022571817]

plot(t,greedy_arr, 'LineWidth', 2)




hold on

fantom_arr = [48.12887666211999,49.79459398377656,50.52411362190529,54.08164106920229,58.240619655153004,61.383778610304816,61.59352937040007,63.87672943247916,65.8976921036234,68.64203480213021]

plot(t,fantom_arr, 'LineWidth', 2)




hold on

sproutpp_arr = [49.1673734019251495,53.2341247494827206,56.4285643303768212,58.191942014152405,60.8824191659898984,61.2215953615274858,63.621603229910265,69.573127540831164,73.912846467880525,76.088277088554166]

plot(t,sproutpp_arr, 'LineWidth', 2)




axis([0.64 1.0 40 80])




hold on



% legend('DSSGS', 'RP\_Greedy', 'Greedy', 'FANTOM', 'SPROUT++', 'Location', 'NorthWest')

legend('DSSGS', 'RP\_Greedy', 'Greedy', 'FANTOM', 'SPROUT++', 'Location', 'SouthEast')

set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)


set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

set(gca,'XTick',t)

ylabel('Objective value','FontName', 'Times New Roman', 'FontSize', 26)


xlabel('Knapsack budget','FontName', 'Times New Roman', 'FontSize', 20)

% set(gca,'Position',[.12 .17 .85 .8])


box on










