classdef InteractiveUI < matlab.apps.AppBase
    % Simple way to create a UI with multiple sliders and 
    % a callback to respond to the update of any slider
    % 
    % Example:
    %
    % figure
    % ax = gca;
    % 
    % s = SliderUI('geometry3D',@(app) (fun(app,ax)));
    % s.AddControl('freq',1, [1 10])
    % s.AddControl('phase',0, [0 360])
    % s.AddControl('amplitude',1, [0 10])
    % s.Open();
    % 
    % function fun(app, ax)
    %   cla(ax)
    %   x = 0:0.001:1;
    %   plot(sin(x*2*pi*app.Values.freq+deg2rad(app.Values.phase))*app.values.amplitude);
    % end


    % Properties that correspond to app components
    properties (Access = public)
        UIFigure        matlab.ui.Figure
        GridLayout      matlab.ui.container.GridLayout
        FigureMenu
        HelpMenu

        % timer
        t

        Values = struct();
        InitialValues
        period = 0.1;

        Data
        HelpText = ''
    end
    properties(Access = private)
        sliderCount = 0;
        updateCallback;
        updating = 0;

        rowHeight = 50;
        rowSpacing = 0;
    end

    % Component initialization
    methods (Access = private)

        function UpdateTimer(app)
            if ( app.updating)
               return
            end
            try
                app.updating = 1;

                if ( app.UIFigure.Visible == "on")
                    app.updateCallback(app);

                    sliders = fieldnames(app.Values);
                    for i=1:length(sliders)
                        slider = sliders{i};

                        if (length(app.GridLayout.Children(i).Children) > 3)
                            % if the control is a slider we need to keep
                            % labels and sliders consistent
                            sliderNum = 3;
                            numFieldNum = 4;

                            if (  class(app.GridLayout.Children(i).Children(sliderNum)) == "matlab.ui.control.Slider" )
                                lims = app.GridLayout.Children(i).Children(sliderNum).Limits;
                                value = app.Values.(slider);
                                value = max(min(value, lims(2)),lims(1));

                                app.GridLayout.Children(i).Children(sliderNum).Value = value;
                                app.GridLayout.Children(i).Children(numFieldNum).Value = value;
                            end
                        end
                    end

                end
                app.updating = 0;
            catch ex
                ex.getReport()
            end
        end

        function Update(app)

        end

        % Value changing function: Slider
        function SliderValueChanging(app, event)
            changingValue = event.Value;
            app.Values.(event.Source.Tag) = changingValue;

            app.Update();
        end

        % Value changed function: Slider
        function SliderValueChanged(app, event)
            value = event.Source.Value;
            app.Values.(event.Source.Tag) = value;
            app.Update();
        end

        % Button pushed function: Button_2
        function ButtonPushed(app, event)
            switch(event.Source.Text)
                case '+'
                    app.Values.(event.Source.Tag) = app.Values.(event.Source.Tag) + 1;
                case '-'
                    app.Values.(event.Source.Tag) = app.Values.(event.Source.Tag) - 1;
            end
            app.Update();
        end

        % Value changed function: EditField
        function EditFieldValueChanged(app, event)
            value = event.Source.Value;
            app.Values.(event.Source.Tag) = value;
            app.Update();
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = InteractiveUI(name, newUpdateCallback, period, helpText)

            if (exist('helpText','var'))
                app.HelpText = helpText;
            end

            % Create UIFigure and components
            N = 2;

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 645 ((app.rowHeight+app.rowSpacing)*N)];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Resize = 'off';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = app.rowHeight;
            app.GridLayout.RowSpacing = app.rowSpacing;

            app.FigureMenu = uimenu(app.UIFigure, 'Text','Menu');
            mitem = uimenu(app.FigureMenu,'Text','Reset');
            mitem.MenuSelectedFcn = @(src,event) app.Reset();

            app.HelpMenu = uimenu(app.UIFigure, 'Text','Help');
            mitem = uimenu(app.HelpMenu,'Text','Show help');
            mitem.MenuSelectedFcn = @(src,event) msgbox(app.HelpText);

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end

            if ( exist("period","var"))
                app.period = period;
            end
            app.updateCallback = newUpdateCallback;
            app.UIFigure.Name = name;

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
            app.InitialValues = app.Values;

            % Setup timer to update psychtoolbox window
            % this is instead of a typical while loop.
            app.t = timer;
            app.t.TimerFcn = @(~,thisevent)app.UpdateTimer;
            app.t.Period = app.period;
            app.t.ExecutionMode = 'fixedRate';
            app.t.TasksToExecute = 100000;

            start(app.t)
        end

        function Reset(app)
            app.Values = app.InitialValues;
            app.Update();
        end

        function Panel = AddPanel(app)
            if (app.sliderCount >= 20)
                error('NO MORE CONTROLS ALLLOWED (20 max)')
            end
            
            app.sliderCount =  app.sliderCount+1;

            % update figure size
            app.UIFigure.Position = [app.UIFigure.Position(1:3) ((app.rowHeight+app.rowSpacing)*app.sliderCount)];

            % Create Panel
            Panel = uipanel(app.GridLayout);
            Panel.Layout.Row = app.sliderCount;
            Panel.Layout.Column = 1;
            % Panel.BorderType = 'none';
        end

        function AddSlider(app, name, defaultvalue, range)

            Panel = app.AddPanel();

            textname = name;
            name = matlab.lang.makeValidName(name);
            app.Values.(name) = defaultvalue;

            % Create EditFieldLabel
            EditFieldLabel = uilabel(Panel);
            EditFieldLabel.HorizontalAlignment = 'right';
            EditFieldLabel.Position = [2 16 140 22];
            EditFieldLabel.Text = textname;

            % Create EditField
            EditField = uieditfield(Panel, 'numeric');
            EditField.Position = [145 16 50 22];
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


        function AddDropDown(app, name, defaultvalue, values)

            Panel = app.AddPanel();

            textname = name;
            name = matlab.lang.makeValidName(name);
            app.Values.(name) = values(defaultvalue);

            % Create EditFieldLabel
            EditFieldLabel = uilabel(Panel);
            EditFieldLabel.HorizontalAlignment = 'right';
            EditFieldLabel.Position = [8 16 115 22];
            EditFieldLabel.Text = textname;

            % Create DropDown
            DropDown = uidropdown(Panel,"Items",values);
            DropDown.Position = [250 14 303 33];
            DropDown.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            DropDown.Tag = name;
            DropDown.Value = DropDown.Items{defaultvalue};

            app.Update();

        end

        function AddMenu(app, name, callback)

            mitem = uimenu(app.FigureMenu,'Text',name);
            mitem.MenuSelectedFcn = @(src,event)MenuUpdate(app, callback);
        end

        function MenuUpdate(app, callback)

            callback(app);
            app.Update();
        end
    end
end