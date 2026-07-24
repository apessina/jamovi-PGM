
module.exports = {

    view_loaded: function(ui, event) {

        makeButtonCompact(ui.duplicateSetup);
        makeButtonCompact(ui.removeSetup);

        ui.duplicateSetup.$el
            .off('click.gcurve')
            .on('click.gcurve', () => {
                this.duplicateSetup(ui);
            });

        ui.removeSetup.$el
            .off('click.gcurve')
            .on('click.gcurve', () => {
                this.removeSetup(ui);
            });
    },


    deps_changed: function(ui, event) {

        const deps = asArray(ui.deps.value());
        const current = asArray(ui.modeling.value());

        const validVariables = new Set(deps);

        const next = current.filter(setup =>
            validVariables.has(setup.var)
        );

        for (const varName of deps) {

            const hasSetup = next.some(setup =>
                setup.var === varName
            );

            if (!hasSetup) {
                next.push(
                    createDefaultSetup(varName, next)
                );
            }
        }

        if (!sameModeling(current, next))
            ui.modeling.setValue(next);
    },


    duplicateSetup: function(ui) {

        const current = asArray(ui.modeling.value());
        const selectedIndex =
            getSelectedModelingIndex(ui);

        if (selectedIndex === -1) {
            return;
        }

        const selected = current[selectedIndex];

        if (selected === undefined)
            return;

        const duplicate = {
            ...selected,
            setupId: makeSetupId(),
            setupLabel: nextSetupLabel(
                selected.var,
                current
            )
        };

        const next = [
            ...current.slice(0, selectedIndex + 1),
            duplicate,
            ...current.slice(selectedIndex + 1)
        ];

        ui.modeling.setValue(next);

        setTimeout(() => {
            ui.modeling.setSelectedRowIndices([
                selectedIndex + 1
            ]);
        }, 0);
    },


    removeSetup: function(ui) {

        const current = asArray(ui.modeling.value());
        const selectedIndex =
            getSelectedModelingIndex(ui);

        if (selectedIndex === -1) {
            return;
        }

        const selected = current[selectedIndex];

        if (selected === undefined)
            return;

        const setupCount = current.filter(setup =>
            setup.var === selected.var
        ).length;

        if (setupCount <= 1) {
            return;
        }

        const next = current.filter(
            (_, index) => index !== selectedIndex
        );

        ui.modeling.setValue(next);

        setTimeout(() => {

            if (next.length === 0)
                return;

            const nextIndex = Math.min(
                selectedIndex,
                next.length - 1
            );

            ui.modeling.setSelectedRowIndices([
                nextIndex
            ]);
        }, 0);
    }

};

// Helpers

function asArray(value) {

    return Array.isArray(value)
        ? value.slice()
        : [];
}


function createDefaultSetup(varName, current) {

    return {
        setupId: makeSetupId(),
        setupLabel: nextSetupLabel(
            varName,
            current
        ),
        var: varName,
        model: 'richards',
        errorVarType: 'none',
        errorCorType: 'none'
    };
}


function makeSetupId() {

    const randomPart = Math.random()
        .toString(36)
        .slice(2, 9);

    return `setup_${Date.now()}_${randomPart}`;
}


function nextSetupLabel(varName, current) {

    const numbers = current
        .filter(setup => setup.var === varName)
        .map(setup => {

            const match = /^M(\d+)$/i.exec(
                setup.setupLabel || ''
            );

            return match === null
                ? 0
                : Number(match[1]);
        });

    const largest = numbers.length === 0
        ? 0
        : Math.max(...numbers);

    return `M${largest + 1}`;
}


function sameModeling(a, b) {

    if (a.length !== b.length)
        return false;

    return a.every((left, index) => {

        const right = b[index];

        return (
            left.setupId === right.setupId &&
            left.setupLabel === right.setupLabel &&
            left.var === right.var &&
            left.model === right.model &&
            left.errorVarType === right.errorVarType &&
            left.errorCorType === right.errorCorType
        );
    });
}


function getSelectedModelingIndex(ui) {

    const selected =
        ui.modeling.getSelectedRowIndices();

    if (
        !Array.isArray(selected) ||
        selected.length === 0
    ) {
        return -1;
    }

    return selected[0];
}


function makeButtonCompact(control) {

    if (!control || !control.$el)
        return;

    control.$el.css({
        width: 'auto',
        minWidth: '0',
        height: '26px',
        padding: '2px 9px',
        fontSize: '12px',
        lineHeight: '18px',
        flexGrow: '0'
    });
}

