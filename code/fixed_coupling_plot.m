V(:, :, :) = r.resultArray(1:meshing, 1:meshing, 1:meshing, chosenCouplingIndex);

for i=1:length(isoContours)
    p = patch(isosurface(taus, pulses, deltas, V, isoContours(i)));
    set(p, 'FaceColor',colors{i}, 'EdgeColor','none')
    view(-69.0636, 25.2000)
    camlight('headlight')
    hold on
end

h=title(strcat('Fidelity Level Surface for Coupling= ', string(couplings(chosenCouplingIndex))));
set(h,'interpreter','latex','fontsize', 20);
h = xlabel('Tau Spacing');
set(h,'interpreter','latex','fontsize',17);
h = ylabel('Pulse Length');
set(h,'interpreter','latex','fontsize',17);
h = zlabel('Chemical Shift');
set(h,'interpreter','latex','fontsize',17);
legend(string(isoContours))
hold off