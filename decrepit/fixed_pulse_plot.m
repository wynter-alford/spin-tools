V(:, :, :) = r.resultArray(chosenPulseIndex, 1:meshing, 1:meshing, 1:meshing);

for i=1:length(isoContours)
    p = patch(isosurface(deltas, taus, couplings, V, isoContours(i)));
    set(p, 'FaceColor',colors{i}, 'EdgeColor','none')
    camlight('right')
    view(-37.5, 30)
    hold on
end
h=title(strcat('Fidelity Level Surface for Pulse= ', string(pulses(chosenPulseIndex))));
set(h,'interpreter','latex','fontsize', 20);
h = xlabel('Chemical Shift');
set(h,'interpreter','latex','fontsize',17);
h = ylabel('Tau Spacing');
set(h,'interpreter','latex','fontsize',17);
h = zlabel('Coupling Strength');
set(h,'interpreter','latex','fontsize',17);
legend(string(isoContours))
hold off
