class SimpleInterface:

    def __init__(self) -> None:
        """ This class is a simple interface to call external programs.

        This class is a simple interface to call external programs. It is designed to be subclassed
        and the method attribute set to the desired program. The call method will then call the program
        with the specified arguments.
        
        """
        return
    
    def set_method(self, method):
        self.method = method
        return
    
    def call(self, **kwargs):
        import subprocess
        dry_run = False
        if "dry_run" in kwargs:
            dry_run = kwargs['dry_run']
            del kwargs['dry_run']
    
        command = [self.method]
        shell=False
        for key, value in kwargs.items():
            if key is "inp_pipe":
                command.extend([f'<', str(value)])
                shell=True
            elif key is "out_pipe":
                command.extend([f'>', str(value)])
                shell=True
            else:
                if value is not None:
                    command.extend([f'-{key}', str(value)])

        if dry_run:
            print(command)
        else:
            print("Executing command")
            proc = subprocess.run(command, shell=shell, encoding='utf-8', stdout=subprocess.PIPE)
            for line in proc.stdout.split('\n'):
                print(f"[{command[0]}] -> {line}")
            print(f"Command Executed:")
            print(f"{' '.join(command)}")
        return
    
class AddToBox(SimpleInterface):
    def __init__(self) -> None:
        """ This class is an interface to the 'add_to_box' program.

        This class is an interface to the 'add_to_box' program. It is a subclass of SimpleInterface
        and the method attribute is set to 'add_to_box'.
        
        """
        self.set_method("AddToBox")
        return
    
    
class Leap(SimpleInterface):
    def __init__(self) -> None:
        """ This class is a simple interface to call the Leap program. """
        self.set_method('tleap')
        return