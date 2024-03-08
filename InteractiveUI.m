classdef InteractiveUI < matlab.apps.AppBase
    % Simple way to create a UI with multiple sliders and 
    % a callback to respond to the update of any slider
    % 
    % Example:
    %
    % figure
    % ax = gca;
    % 
    % s = SliderUI('geometry3D',@(sliderValues) (fun(sliderValues,ax)));
    % s.AddControl('freq',1, [1 10])
    % s.AddControl('phase',0, [0 360])
    % s.AddControl('amplitude',1, [0 10])
    % s.Open();
    % 
    % function fun(sliderValues, ax)
    %   cla(ax)
    %   x = 0:0.001:1;
    %   plot(sin(x*2*pi*values.freq+deg2rad(sliderValues.phase))*values.amplitude);
    % end


    % Properties that correspond to app components
    properties (Access = public)
        UIFigure        matlab.ui.Figure
        GridLayout      matlab.ui.container.GridLayout

        % timer
        t

        SliderValues = struct();
    end
    properties(Access = private)
        sliderCount = 0;
        updateCallback;
        updating = 0;
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 645 721];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Resize = 'off';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = {60, 60, 60, 60, 60, 60, 60, 60, 60, 60};
        end

        function InitTimer(app, period)
            if ( ~exist('period','var'))
                period = 0.1;
            end

            % Setup timer to update psychtoolbox window
            % this is instead of a typical while loop.
            app.t = timer;
            app.t.TimerFcn = @(~,thisevent)app.UpdateTimer;
            app.t.Period = period;
            app.t.ExecutionMode = 'fixedRate';
            app.t.TasksToExecute = 100000;

            start(app.t)
        end

        function UpdateTimer(app)
            if ( app.updating)
               return
            end
            app.updating = 1;

            if ( app.UIFigure.Visible == "on")
                app.updateCallback(app.SliderValues);

                sliders = fieldnames(app.SliderValues);
                for i=1:length(sliders)
                    slider = sliders{i};

                    sliderNum = 3;
                    numFieldNum = 4;

                    lims = app.GridLayout.Children(i).Children(sliderNum).Limits;
                    value = app.SliderValues.(slider);
                    value = max(min(value, lims(2)),lims(1));

                    app.GridLayout.Children(i).Children(sliderNum).Value = value;
                    app.GridLayout.Children(i).Children(numFieldNum).Value = value;
                end

            end
            app.updating = 0;
        end

        function Update(app)

        end

        % Value changing function: Slider
        function SliderValueChanging(app, event)
            changingValue = event.Value;
            app.SliderValues.(event.Source.Tag) = changingValue;

            app.Update();
        end

        % Value changed function: Slider
        function SliderValueChanged(app, event)
            value = event.Source.Value;
            app.SliderValues.(event.Source.Tag) = value;
            app.Update();
        end

        % Button pushed function: Button_2
        function ButtonPushed(app, event)
            switch(event.Source.Text)
                case '+'
                    app.SliderValues.(event.Source.Tag) = app.SliderValues.(event.Source.Tag) + 1;
                case '-'
                    app.SliderValues.(event.Source.Tag) = app.SliderValues.(event.Source.Tag) - 1;
            end
            app.Update();
        end

        % Value changed function: EditField
        function EditFieldValueChanged(app, event)
            value = event.Source.Value;
            app.SliderValues.(event.Source.Tag) = value;
            app.Update();
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = InteractiveUI(name, newUpdateCallback, period)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end

            app.updateCallback = newUpdateCallback;
            app.UIFigure.Name = name;

            app.InitTimer(period);
        end

        % Code that executes before app deletion
        function delete(app)

            stop(app.t)
            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end

        function Open(app)

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end

        function AddControl(app, name, defaultvalue, range)

            if (app.sliderCount >= 10)
                error('NO MORE SLIDERS ALLLOWED')
            end
            app.sliderCount =  app.sliderCount+1;

            textname = name;
            name = matlab.lang.makeValidName(name);
            app.SliderValues.(name) = defaultvalue;

            % Create Panel
            Panel = uipanel(app.GridLayout);
            Panel.Layout.Row = app.sliderCount;
            Panel.Layout.Column = 1;

            % Create EditFieldLabel
            EditFieldLabel = uilabel(Panel);
            EditFieldLabel.HorizontalAlignment = 'right';
            EditFieldLabel.Position = [8 16 115 22];
            EditFieldLabel.Text = textname;

            % Create EditField
            EditField = uieditfield(Panel, 'numeric');
            EditField.Position = [133 16 60 22];
            EditField.ValueChangedFcn = createCallbackFcn(app, @EditFieldValueChanged, true);
            EditField.Tag = name;

            % Create Slider
            Slider = uislider(Panel);
            Slider.Position = [250 34 303 3];
            Slider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            Slider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanging, true);
            Slider.Tag = name;
            Slider.Limits = range;

            % Create Button
            Button = uibutton(Panel, 'push');
            Button.Position = [574 14 31 23];
            Button.Text = '+';
            Button.ButtonPushedFcn = createCallbackFcn(app, @ButtonPushed, true);
            Button.Tag = name;

            % Create Button_2
            Button_2 = uibutton(Panel, 'push');
            Button_2.ButtonPushedFcn = createCallbackFcn(app, @ButtonPushed, true);
            Button_2.Position = [204 16 29 23];
            Button_2.Text = '-';
            Button_2.Tag = name;

            app.Update();

        end
    end
end