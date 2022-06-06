import click

class SectionedHelpGroup(click.Group):
    """Sections commands into help groups"""

    def __init__(self, *args, **kwargs):
        self.grouped_commands = kwargs.pop('grouped_commands', {})
        commands = {}
        for group, command_list in self.grouped_commands.items():
            for cmd in command_list:
                cmd.help_group = group
                commands[cmd.name] = cmd

        super(SectionedHelpGroup, self).__init__(
            *args, commands=commands, **kwargs)

    def command(self, *args, **kwargs):
        help_group = kwargs.pop('help_group')
        decorator = super(SectionedHelpGroup, self).command(*args, **kwargs)

        def new_decorator(f):
            cmd = decorator(f)
            cmd.help_group = help_group
            self.grouped_commands.setdefault(help_group, []).append(cmd)
            return cmd

        return new_decorator

    def format_commands(self, ctx, formatter):
        for group, cmds in self.grouped_commands.items():
            rows = []
            for subcommand in self.list_commands(ctx):
                cmd = self.get_command(ctx, subcommand)
                if cmd is None or cmd.help_group != group:
                    continue
                help = cmd.get_help(ctx).split('\n\n')[1].strip() or ''
                if len(help) > 57:
                    help = help[:57].strip() + "..."
                rows.append((subcommand, help))

            if rows:
                with formatter.section(group):
                    formatter.write_dl(rows)
