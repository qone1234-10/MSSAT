function plot_age_curves_colormap(ax, years, ages, Xay, cmapAge, ttl, isJoint, emphBand)
    axes(ax); %#ok<LAXES>
    hold(ax,'on');

    ax.Box = 'on';
    ax.LineWidth = 1.0;
    ax.FontName = 'Arial';
    ax.FontSize = 11;
    ax.XGrid = 'on'; ax.YGrid = 'on';
    ax.GridAlpha = 0.15;
    ax.TickDir = 'out';

    emph_idx = find(ages>=emphBand(1) & ages<=emphBand(2));

    lw_base = 0.9;
    lw_emph = 1.8;

    for k = 1:numel(ages)
        lw = lw_base;
        if any(k == emph_idx)
            lw = lw_emph;
        end
        plot(ax, years, Xay(k,:), 'Color', cmapAge(k,:), 'LineWidth', lw);
    end

    xlabel(ax,'Year');
    ylabel(ax,'Value');
    title(ax, ttl, 'FontWeight','bold');
    xlim(ax, [min(years) max(years)]);

    if isJoint
        yl = quantile(Xay(:), [0.01 0.99]);
        pad = 0.15 * range(yl);
    else
        yl = quantile(Xay(:), [0.02 0.98]);
        pad = 0.25 * range(yl);
    end
    ylim(ax, [yl(1)-pad, yl(2)+pad]);

    hold(ax,'off');
end