
module.exports = {

    view_loaded: function(ui, event) {
        syncParameterRows(ui);
    },
    
    view_updated: function(ui, event) {
        syncParameterRows(ui);
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
    
    duplicateSetup_clicked: function(ui, event) {

        if (ui.duplicateSetup.value() !== true)
            return;

        this.duplicateSetup(ui);
    },

    removeSetup_clicked: function(ui, event) {

        if (ui.removeSetup.value() !== true)
            return;

        this.removeSetup(ui);
    },

    modelParam_changed: function(ui, event) {
        syncParameterRows(ui);
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
        fixValue: '',
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
            left.errorCorType === right.errorCorType && 
            left.fixValue === right.fixValue
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


function syncParameterRows(ui) {

    const names = [
        'A',
        'K',
        ui.modelParam.value() === 'W0' ? 'W0' : 'Ti',
        'd'
    ];

    const fixedCurrent = asArray(ui.fixedTerms.value());
    const randomCurrent = asArray(ui.randomGroups.value());

    const fixedNext =
        normalizeParameterRows(fixedCurrent, names);

    const randomNext =
        normalizeParameterRows(randomCurrent, names);

    ui.view.model.options.beginEdit();

    if (!sameValue(fixedCurrent, fixedNext))
        ui.fixedTerms.setValue(fixedNext);

    if (!sameValue(randomCurrent, randomNext))
        ui.randomGroups.setValue(randomNext);

    ui.view.model.options.endEdit();
}


function normalizeParameterRows(current, names) {

    return names.map((param, index) => {

        const row = current[index];

        if (
            row &&
            !Array.isArray(row) &&
            typeof row === 'object'
        ) {
            return {
                param: param,
                terms: asArray(row.terms)
            };
        }

        return {
            param: param,
            terms: Array.isArray(row)
                ? row.slice()
                : []
        };
    });
}


function sameValue(a, b) {
    return JSON.stringify(a) === JSON.stringify(b);
}

